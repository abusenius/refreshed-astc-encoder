#ifndef ASTC_OPENCL_INCLUDED

#define ASTC_OPENCL_INCLUDED

#include <CL/opencl.h>
#include <assert.h>
#include <stdio.h>

extern cl_device_id opencl_device;
extern cl_program opencl_program;
extern cl_context opencl_context;

#define OPENCL_KERNELS_SOURCE_PATH ".\\"
#define OPENCL_KERNEL_FILES "astc_find_best_partitioning.cl", "astc_averages_and_directions.cl"
#define OPENCL_COMPILER_OPTIONS "-cl-mad-enable"

#define OCL_IS_BLOCKING false

#define OCL_CHECK_STATUS(str) if(status != CL_SUCCESS) { fprintf(stderr, "%s, errorcode: %i %s\n", str, status, cl_errcode_to_str(status)); exit(-1); }
#define OCL_RELEASE_OBJECT(type, name) { status = clRelease##type(name); OCL_CHECK_STATUS("Cannot release "#type" "#name); }
#define OCL_CREATE_BUFFER(name, flags, size, source) { name = clCreateBuffer(opencl_context, flags, size, source, &status);\
								OCL_CHECK_STATUS("Cannot create buffer "#name); }
#define OCL_WRITE_BUFFER(name, size, src) { status = clEnqueueWriteBuffer(opencl_queue, name, OCL_IS_BLOCKING, 0, size, src, 0, NULL, NULL);\
								OCL_CHECK_STATUS("Cannot copy data to bufer "#name); }
#define OCL_READ_BUFFER(name, size, dest) { status = clEnqueueReadBuffer(opencl_queue, name, OCL_IS_BLOCKING, 0, size, dest, 0, NULL, NULL);\
								OCL_CHECK_STATUS("Cannot read data from bufer "#name); }
#define OCL_MAP_BUFFER(name, ptr_type, ptr, flags, size) { ptr = static_cast<ptr_type> (clEnqueueMapBuffer(opencl_queue, name, OCL_IS_BLOCKING, flags, 0, size, 0, NULL, NULL, &status)); \
								OCL_CHECK_STATUS("Error in clEnqueueMapBuffer "#name); }
#define OCL_UNMAP_BUFFER(name, ptr) { status = clEnqueueUnmapMemObject(opencl_queue, name, ptr, 0, NULL, NULL);\
								OCL_CHECK_STATUS("Error in clEnqueueUnmapMemObject "#name); }


#define OCL_CREATE_KERNEL(module, name) { module.##name = clCreateKernel(opencl_program, #name, &status);\
								OCL_CHECK_STATUS("Cannot create kernel "#name); }
#define OCL_SET_KERNEL_ARG(kernname, seq, argname) { status = set_kernel_argument(kernname, seq, argname);\
								OCL_CHECK_STATUS("Cannot set argument "#seq" ("#argname") for kernel "#kernname); }



void init_opencl(cl_uint platform_number, cl_uint device_number, int silentmode, int batch_size, int xdim, int ydim, int zdim, const error_weighting_params* ewp, astc_decode_mode decode_mode);
void destroy_opencl();
char const* cl_errcode_to_str(cl_int status);

enum class ocl_buffer_type
{
	DEVICE,
	DEVICE_PREPINNED,
	PINNED_PAIR
};

// memory objects used by both host and OpenCL device
template <typename T, ocl_buffer_type BUF_TYPE> class ocl_buffer
{
public:
	ocl_buffer() : valid(false) {};
	~ocl_buffer();

	cl_mem dev_buf;
	cl_mem host_buf;
	T * host_ptr;

	void create_buffer(cl_mem_flags mem_flags, size_t count, cl_command_queue _opencl_queue);
	void write_to_device();
	void read_from_device();
	T& operator[](size_t idx) { return host_ptr[idx]; };
	const T& operator[](size_t idx) const { return host_ptr[idx]; };

private:
	size_t size;
	bool valid;
	cl_command_queue opencl_queue;
};


template<typename T, ocl_buffer_type BUF_TYPE>
inline void ocl_buffer<T, BUF_TYPE>::create_buffer(cl_mem_flags mem_flags, size_t count, cl_command_queue _opencl_queue)
{
	cl_int status;

	size = sizeof(T) * count;
	opencl_queue = _opencl_queue;

	switch (BUF_TYPE)
	{
	case ocl_buffer_type::DEVICE:
		host_ptr = new T[count];
		OCL_CREATE_BUFFER(dev_buf, mem_flags, size, NULL);
		break;
	case ocl_buffer_type::DEVICE_PREPINNED:
		OCL_CREATE_BUFFER(dev_buf, mem_flags | CL_MEM_ALLOC_HOST_PTR, size, NULL);
		OCL_MAP_BUFFER(dev_buf, T*, host_ptr, CL_MAP_WRITE_INVALIDATE_REGION, size);
		break;
	case ocl_buffer_type::PINNED_PAIR:
		OCL_CREATE_BUFFER(dev_buf, mem_flags | CL_MEM_HOST_NO_ACCESS, size, NULL);
		OCL_CREATE_BUFFER(host_buf, mem_flags | CL_MEM_ALLOC_HOST_PTR, size, NULL);
		OCL_MAP_BUFFER(host_buf, T*, host_ptr, CL_MAP_WRITE_INVALIDATE_REGION, size);
		break;
	default:
		ASTC_CODEC_INTERNAL_ERROR;
		break;
	}

	status = clRetainCommandQueue(opencl_queue);
	assert(status == CL_SUCCESS);
	valid = true;
}


template<typename T, ocl_buffer_type BUF_TYPE>
inline void ocl_buffer<T, BUF_TYPE>::write_to_device()
{
	cl_int status;

	switch (BUF_TYPE)
	{
	case ocl_buffer_type::DEVICE:
		OCL_WRITE_BUFFER(dev_buf, size, host_ptr);
		break;
	case ocl_buffer_type::DEVICE_PREPINNED:
		OCL_UNMAP_BUFFER(dev_buf, host_ptr);
		break;
	case ocl_buffer_type::PINNED_PAIR:
		OCL_WRITE_BUFFER(dev_buf, size, host_ptr);
		break;
	default:
		ASTC_CODEC_INTERNAL_ERROR;
		break;
	}
}

template<typename T, ocl_buffer_type BUF_TYPE>
inline void ocl_buffer<T, BUF_TYPE>::read_from_device()
{
	cl_int status;

	switch (BUF_TYPE)
	{
	case ocl_buffer_type::DEVICE:
		OCL_READ_BUFFER(dev_buf, size, host_ptr);
		break;
	case ocl_buffer_type::DEVICE_PREPINNED:
		OCL_MAP_BUFFER(dev_buf, T*, host_ptr, CL_MAP_READ | CL_MAP_WRITE, size);
		break;
	case ocl_buffer_type::PINNED_PAIR:
		OCL_WRITE_BUFFER(dev_buf, size, host_ptr);
		break;
	default:
		ASTC_CODEC_INTERNAL_ERROR;
		break;
	}
}


template<typename T, ocl_buffer_type BUF_TYPE>
ocl_buffer<T, BUF_TYPE>::~ocl_buffer()
{
	if (!valid)
		return;

	cl_int status;

	switch (BUF_TYPE)
	{
	case ocl_buffer_type::DEVICE:
		OCL_RELEASE_OBJECT(MemObject, dev_buf);
		delete[] host_ptr;
		break;
	case ocl_buffer_type::DEVICE_PREPINNED:
		OCL_UNMAP_BUFFER(dev_buf, host_ptr);
		OCL_RELEASE_OBJECT(MemObject, dev_buf);
		break;
	case ocl_buffer_type::PINNED_PAIR:
		OCL_UNMAP_BUFFER(host_buf, host_ptr);
		OCL_RELEASE_OBJECT(MemObject, host_buf);
		OCL_RELEASE_OBJECT(MemObject, dev_buf);
		break;
	default:
		ASTC_CODEC_INTERNAL_ERROR;
		break;
	}

	OCL_RELEASE_OBJECT(CommandQueue, opencl_queue);

	valid = false;
}


template <typename T, ocl_buffer_type BT> inline cl_int set_kernel_argument(cl_kernel kernel, cl_uint seq, const ocl_buffer<T, BT>& arg)
{
	return clSetKernelArg(kernel, seq, sizeof(cl_mem), &arg.dev_buf);
}

template <typename T> inline cl_int set_kernel_argument(cl_kernel kernel, cl_uint seq, const T& arg)
{
	return clSetKernelArg(kernel, seq, sizeof(T), &arg);
}


#endif
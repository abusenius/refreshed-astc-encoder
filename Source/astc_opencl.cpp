#include "astc_codec_batch.h"
#include "metrohash64.h"


#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <string>
#include <vector>


cl_platform_id opencl_platform;
cl_device_id opencl_device;
cl_program opencl_program;
cl_context opencl_context;

#ifdef WIN32
	#define PATH_SEPPARATOR "\\"
#else
	#define PATH_SEPPARATOR "/"
#endif // WIN32

static cl_int printDeviceInfo(cl_device_id device)
{
	cl_int status = CL_SUCCESS;
	char buf[512], buf2[512];
	cl_uint clui, clui2;
	cl_ulong clul, clul2;
	cl_device_local_mem_type memtype;
	char *memtype_str[4] = {"none", "local", "global", "unknown"};
	size_t  clsz;

	status = clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(buf), buf, NULL);
	if (status == CL_SUCCESS) printf("\t\tVendor: %s\n", buf);

	status = clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(buf), buf, NULL);
	if (status == CL_SUCCESS) printf("\t\tName: %s\n", buf);

	status = clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(buf), buf, NULL);
	if (status == CL_SUCCESS) printf("\t\tDriver Version: %s\n", buf);

	status = clGetDeviceInfo(device, CL_DEVICE_VERSION, sizeof(buf), buf, NULL);
	status = clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, sizeof(buf2), buf2, NULL);
	if (status == CL_SUCCESS) printf("\t\tDev/OpenCL C ver: %s/%s\n", buf, buf2);

	status = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &clui, NULL);
	status |= clGetDeviceInfo(device, CL_DEVICE_PARTITION_MAX_SUB_DEVICES, sizeof(cl_uint), &clui2, NULL);
	if (status == CL_SUCCESS) printf("\t\tUnits/MaxPartitions: %u / %u\n", clui, clui2);

	status = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &clui, NULL);
	if (status == CL_SUCCESS) printf("\t\tMax dim: %u\n", clui);

	status = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &clsz, NULL);
	if (status == CL_SUCCESS) printf("\t\tWG size: %u\n", (cl_uint)clsz);

	status = clGetDeviceInfo(device, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint), &clui, NULL);
	if (status == CL_SUCCESS) printf("\t\tFreqs: %u\n", clui);

	status = clGetDeviceInfo(device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), &clul, NULL);
	status |= clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &clul2, NULL);
	if (status == CL_SUCCESS) printf("\t\tGlobal Mem size/max alloc: %u MiB/ %u MiB\n", (cl_uint)(clul2>>20), (cl_uint)(clul >> 20));

	status = clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(cl_device_local_mem_type), &memtype, NULL);
	status |= clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &clul, NULL);
	if (status == CL_SUCCESS) printf("\t\tLocal Mem size(type): %u KiB (%s)\n", (cl_uint)(clul >> 10), memtype_str[MIN(memtype, 3)]);

	status = clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, sizeof(cl_uint), &clui, NULL);
	status |= clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, sizeof(cl_ulong), &clul, NULL);
	if (status == CL_SUCCESS) printf("\t\tGlobal cache size/line size: %u KiB / %u B\n", (cl_uint)(clul >> 10), clui);

	return CL_SUCCESS;
}

static cl_int printPlatformInfo(cl_platform_id platform)
{
	cl_int status = CL_SUCCESS;
	char buf[512];

	status = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, sizeof(buf), buf, NULL);
	if (status == CL_SUCCESS) printf("\tVendor: %s\n", buf);

	status = clGetPlatformInfo(platform, CL_PLATFORM_NAME, sizeof(buf), buf, NULL);
	if (status == CL_SUCCESS) printf("\tName: %s\n", buf);

	status = clGetPlatformInfo(platform, CL_PLATFORM_PROFILE, sizeof(buf), buf, NULL);
	if (status == CL_SUCCESS) printf("\tProfile: %s\n", buf);

	status = clGetPlatformInfo(platform, CL_PLATFORM_VERSION, sizeof(buf), buf, NULL);
	if (status == CL_SUCCESS) printf("\tVersion: %s\n", buf);

	return CL_SUCCESS;
}

static cl_int readFile(const char *path, const char *filename, char** sourceString, size_t *file_size)
{
	FILE *fp;
	size_t err;
	size_t size;

	std::string fullFilename(path);
	fullFilename += PATH_SEPPARATOR;
	fullFilename += filename;

	char *source;
	fopen_s(&fp, fullFilename.c_str(), "rb");
	if (fp == NULL) {
		fprintf(stderr, "Could not open file (read): %s\n", fullFilename.c_str());
		return CL_INVALID_VALUE;
	}

	err = fseek(fp, 0, SEEK_END);
	if (err != 0) {
		fprintf(stderr, "Error seeking to end of file\n");
		return CL_INVALID_VALUE;
	}

	size = ftell(fp);
	if (size < 0) {
		fprintf(stderr, "Error getting file position\n");
		return CL_INVALID_VALUE;
	}

	err = fseek(fp, 0, SEEK_SET);
	if (err != 0) {
		fprintf(stderr, "Error seeking to start of file\n");
		return CL_INVALID_VALUE;
	}

	source = new char[size + 1];

	err = fread(source, 1, size, fp);
	if (err != size) {
		fprintf(stderr, "only %u bytes read\n", (unsigned int)err);
		delete[] source;
		return CL_INVALID_VALUE;
	}

	source[size] = '\0';

	*sourceString = source;
	if (file_size)
		*file_size = size;
	return CL_SUCCESS;
}



static void generateBinaryFilename(cl_platform_id platform, cl_device_id device, uint64_t source_hash, char *filename)
{
	cl_int status;
	MetroHash64 hasher;
	uint8_t buf[512], plat_name[512];
	size_t size;

	hasher.Initialize();

	status = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, sizeof(buf), buf, &size);
	if (status == CL_SUCCESS) hasher.Update(buf, size);

	status = clGetPlatformInfo(platform, CL_PLATFORM_VERSION, sizeof(buf), buf, &size);
	if (status == CL_SUCCESS) hasher.Update(buf, size);

	status = clGetPlatformInfo(platform, CL_PLATFORM_NAME, sizeof(plat_name), plat_name, &size);
	if (status == CL_SUCCESS) hasher.Update(plat_name, size);

	status = clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(buf), buf, &size);
	if (status == CL_SUCCESS) hasher.Update(buf, size);

	status = clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(buf), buf, &size);
	if (status == CL_SUCCESS) hasher.Update(buf, size);

	status = clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(buf), buf, &size);
	if (status == CL_SUCCESS) hasher.Update(buf, size);

	char * dev_name_short = strtok((char*)buf, " \"*:,<>?|()\\/`'_");
	char * plat_name_short = strtok((char*)plat_name, " \"*:,<>?|()\\/`'_");

	uint64_t driver_hash;
	hasher.Finalize((uint8_t*)&driver_hash);

	sprintf(filename, "%.*s_%.*s_%016llX_%016llX.bin", 20, plat_name_short, 20, dev_name_short, driver_hash, source_hash);
}

static cl_int loadBinary(cl_context context, cl_device_id device, const char *path, const char *fname, cl_program *program, int silentmode)
{
	cl_int status, binary_status;
	char *binary;
	size_t size;
	if (CL_SUCCESS != readFile(path, fname, &binary, &size))
		return CL_INVALID_BINARY;

	if (!silentmode)
		printf("Found precompiled kernels %s\n", fname);

	*program = clCreateProgramWithBinary(context, 1, &device, &size, (const unsigned char**)&binary, &binary_status, &status);
	delete[] binary;

	if (binary_status)
	{
		if (!silentmode)
			printf("clCreateProgramWithBinary() failed, binary_status: %i %s\n", binary_status, cl_errcode_to_str(binary_status));
		return binary_status;
	}

	if (status && !silentmode)
		printf("clCreateProgramWithBinary() failed, errorcode: %i %s\n", status, cl_errcode_to_str(status));

	return status;
}

static cl_int saveBinary(const char *path, const char *fname, cl_program opencl_program, const char *compilerOptions, int silentmode)
{
	cl_int status;


	// retrieve compiled kernel
	size_t size, bytes_written;
	status = clGetProgramInfo(opencl_program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &size, NULL);
	OCL_CHECK_STATUS_R("Unable to retrieve program binary size", -1);
	if (size == 0)
	{
		if (!silentmode)
			printf("The binary is not available for current device");
		return -1;
	}

	std::vector<unsigned char> binary(size);
	unsigned char *ptr = &binary[0];
	status = clGetProgramInfo(opencl_program, CL_PROGRAM_BINARIES, sizeof(void*), &ptr, NULL);
	OCL_CHECK_STATUS_R("Unable to retrieve program binary", -1);
	

	// save binary to file
	std::string fullFilename(path);
	fullFilename += PATH_SEPPARATOR;
	fullFilename += fname;

	FILE *fp = fopen(fullFilename.c_str(), "wb");
	if (!fp)
	{
		fprintf(stderr, "Could not open file (write): %s\n", fullFilename.c_str());
		return -1;
	}

	bytes_written = fwrite(&binary[0], sizeof(char), size, fp);
	fclose(fp);
	if ((size != bytes_written) && !silentmode)
	{
		printf("The binary size is %zi but only %zi bytes written to file\n", size, bytes_written);
		return -1;
	}

	// save info to index.txt
	fullFilename = path;
	fullFilename += PATH_SEPPARATOR;
	fullFilename += "index.txt";

	fp = fopen(fullFilename.c_str(), "a");
	if (!fp)
	{
		fprintf(stderr, "Could not open file (append): %s\n", fullFilename.c_str());
		return -1;
	}

	std::string record(fname);
	record += "\t";
	record += compilerOptions;
	record += "\n";

	size = record.size();
	bytes_written = fwrite(record.c_str(), sizeof(char), size, fp);
	fclose(fp);
	if ((size != bytes_written) && !silentmode)
	{
		printf("The log size is %zi but only %zi bytes written to file\n", size, bytes_written);
		return -1;
	}
	
	return CL_SUCCESS;
}

void init_opencl(const opencl_options * oclo, int batch_size, int xdim, int ydim, int zdim, const error_weighting_params * ewp, astc_decode_mode decode_mode, int silentmode)
{
	cl_uint numPlatforms = 0;
	cl_uint numDevices = 0;
	cl_platform_id *platforms;
	cl_device_id *devices;
	cl_int status;
	MetroHash64 hasher;

	// get platform
	status = clGetPlatformIDs(0, NULL, &numPlatforms);
	OCL_CHECK_STATUS("Cannot get number of platforms");
	if (numPlatforms <= oclo->platform) {
		fprintf(stderr, "Only %i platforms found. Platform with index %i does not exists.", numPlatforms, oclo->platform);
		exit(-1);
	}

	platforms = new cl_platform_id[numPlatforms];

	status = clGetPlatformIDs(numPlatforms, platforms, NULL);
	opencl_platform = platforms[oclo->platform];
	delete[] platforms;

	if (!silentmode) {
		printf("%i platforms found. Using platform %i\n", numPlatforms, oclo->platform);
		printPlatformInfo(opencl_platform);
	}

	// get device
	status = clGetDeviceIDs(opencl_platform, CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
	OCL_CHECK_STATUS("Cannot get number of devices");
	if (numDevices <= oclo->device) {
		fprintf(stderr, "Only %i devices found. Device with index %i does not exists.", numDevices, oclo->device);
		exit(-1);
	}

	devices = new cl_device_id[numDevices];

	status = clGetDeviceIDs(opencl_platform, CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);
	opencl_device = devices[oclo->device];
	delete[] devices;
	if (!silentmode) {
		printf("\n%i devices found on platform %i. Using device %i\n", numDevices, oclo->platform, oclo->device);
		printDeviceInfo(opencl_device);
	}

	// create context
	cl_context_properties context_props[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)opencl_platform, 0};
	opencl_context = clCreateContext(context_props, 1, &opencl_device, NULL, NULL, &status);
	OCL_CHECK_STATUS("Cannot create context");


	// read source files with kernels
	hasher.Initialize();
	char *source_names[] = { OPENCL_KERNEL_FILES };
#define FILE_COUNT ((int)(sizeof(source_names)/sizeof(source_names[0])))
	char *sources[FILE_COUNT];
	for (int i = 0; i < FILE_COUNT; i++)
	{
		size_t filesize;
		status = readFile(oclo->kernels_source_path, source_names[i], &sources[i], &filesize);
		OCL_CHECK_STATUS("Cannot read OpenCL kernel source file");
		hasher.Update((uint8_t*)sources[i], filesize);
	}

	//set kernel compile-time constants
	int texels_per_block = xdim * ydim * zdim;
	float weight_imprecision_estim_squared = calculate_weight_imprecision(texels_per_block);

	char compile_flags[1024] = "";
	if (sizeof(void*) == 4)
		strcat(compile_flags, " -D OCL_USE_32BIT_POINTERS");

	if (sizeof(void*) == 8)
		strcat(compile_flags, " -D OCL_USE_64BIT_POINTERS");

	char compilerOptions[4096];
	setlocale(LC_NUMERIC, "C");
	sprintf(compilerOptions, "%s -I %s -D XDIM=%i -D YDIM=%i -D ZDIM=%i -D TEXELS_PER_BLOCK=%i -D WEIGHT_IMPRECISION_ESTIM_SQUARED=%gf -D PLIMIT=%i"
                              " -D DLIMIT_1PLANE=%i -D DLIMIT_2PLANES=%i -D WLIMIT_1PLANE=%i -D WLIMIT_2PLANES=%i -D BATCH_SIZE=%i %s",
		OPENCL_COMPILER_OPTIONS, ".", xdim, ydim, zdim, texels_per_block, weight_imprecision_estim_squared, ewp->partition_search_limit,
		ewp->decimation_mode_limit_1plane, ewp->decimation_mode_limit_2planes, ewp->weight_mode_limit_1plane, ewp->weight_mode_limit_2planes, batch_size, compile_flags);
	if (!silentmode)
		printf("Batch size: %i\nOpenCL compiler options:\n%s\n\n", batch_size, compilerOptions);

	hasher.Update((const uint8_t*)compilerOptions, strlen(compilerOptions));
	uint64_t source_hash;
	hasher.Finalize((uint8_t*)&source_hash);

	// try to load precompiled binaries from disk
	char fname[512];
	bool loaded_from_file = false;
	generateBinaryFilename(opencl_platform, opencl_device, source_hash, fname);
	if (oclo->enable_cache_read)
	{
		status = loadBinary(opencl_context, opencl_device, oclo->kernels_cache_path, fname, &opencl_program, silentmode);
		loaded_from_file = (CL_SUCCESS == status);
	}
	if (!loaded_from_file)
		opencl_program = clCreateProgramWithSource(opencl_context, FILE_COUNT, (const char **)sources, NULL, &status);
	for (int i = 0; i < FILE_COUNT; i++)
		delete[] sources[i];
	
	OCL_CHECK_STATUS("Cannot create OpenCL program object");
	
	// build program
	if (!silentmode)
		printf(loaded_from_file ? "Building kernels from file... " : "Compiling kernels...  ");
	status = clBuildProgram(opencl_program, 1, &opencl_device, compilerOptions, NULL, NULL);
	if (status != CL_SUCCESS)
	{
		fprintf(stderr, "Cannot build OpenCL program");

		if (!silentmode) {
			//getting log size
			size_t logSize = 1;
			status = clGetProgramBuildInfo(opencl_program, opencl_device, CL_PROGRAM_BUILD_LOG, NULL, NULL, &logSize);
			OCL_CHECK_STATUS("Cannot get build log size");

			//retreiving build log
			std::vector<char> buildLog(logSize);
			
			status = clGetProgramBuildInfo(opencl_program, opencl_device, CL_PROGRAM_BUILD_LOG, logSize, &buildLog[0], NULL);
			if (status != CL_SUCCESS) {
				fprintf(stderr, "Unable to retreive Build log");
				exit(-1);
			}

			fprintf(stderr, "Build Log:\n%.*s\n", (int)logSize, &buildLog[0]);
		}
		getchar();
		exit(-1);
	}

	if (!silentmode)
		printf("Done\n\n");

	if (oclo->enable_cache_write && !loaded_from_file)
		saveBinary(oclo->kernels_cache_path, fname, opencl_program, compilerOptions, silentmode);
}


void destroy_opencl()
{
	cl_int status;

	OCL_RELEASE_OBJECT(Program, opencl_program);
	OCL_RELEASE_OBJECT(Context, opencl_context);
	OCL_RELEASE_OBJECT(Device, opencl_device);
}

char const* cl_errcode_to_str(cl_int status)
{
	static char *errors[] = {
		"CL_SUCCESS",
		"CL_DEVICE_NOT_FOUND",
		"CL_DEVICE_NOT_AVAILABLE",
		"CL_COMPILER_NOT_AVAILABLE",
		"CL_MEM_OBJECT_ALLOCATION_FAILURE",
		"CL_OUT_OF_RESOURCES",
		"CL_OUT_OF_HOST_MEMORY",
		"CL_PROFILING_INFO_NOT_AVAILABLE",
		"CL_MEM_COPY_OVERLAP",
		"CL_IMAGE_FORMAT_MISMATCH",
		"CL_IMAGE_FORMAT_NOT_SUPPORTED",
		"CL_BUILD_PROGRAM_FAILURE",
		"CL_MAP_FAILURE",
		"CL_MISALIGNED_SUB_BUFFER_OFFSET",
		"CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST",
		"CL_COMPILE_PROGRAM_FAILURE",
		"CL_LINKER_NOT_AVAILABLE",
		"CL_LINK_PROGRAM_FAILURE",
		"CL_DEVICE_PARTITION_FAILED",
		"CL_KERNEL_ARG_INFO_NOT_AVAILABLE",
		"CODE_UNKNOWN",
		"CODE_UNKNOWN",
		"CODE_UNKNOWN",
		"CODE_UNKNOWN",
		"CODE_UNKNOWN",
		"CODE_UNKNOWN",
		"CODE_UNKNOWN",
		"CODE_UNKNOWN",
		"CODE_UNKNOWN",
		"CODE_UNKNOWN",
		"CL_INVALID_VALUE",
		"CL_INVALID_DEVICE_TYPE",
		"CL_INVALID_PLATFORM",
		"CL_INVALID_DEVICE",
		"CL_INVALID_CONTEXT",
		"CL_INVALID_QUEUE_PROPERTIES",
		"CL_INVALID_COMMAND_QUEUE",
		"CL_INVALID_HOST_PTR",
		"CL_INVALID_MEM_OBJECT",
		"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
		"CL_INVALID_IMAGE_SIZE",
		"CL_INVALID_SAMPLER",
		"CL_INVALID_BINARY",
		"CL_INVALID_BUILD_OPTIONS",
		"CL_INVALID_PROGRAM",
		"CL_INVALID_PROGRAM_EXECUTABLE",
		"CL_INVALID_KERNEL_NAME",
		"CL_INVALID_KERNEL_DEFINITION",
		"CL_INVALID_KERNEL",
		"CL_INVALID_ARG_INDEX",
		"CL_INVALID_ARG_VALUE",
		"CL_INVALID_ARG_SIZE",
		"CL_INVALID_KERNEL_ARGS",
		"CL_INVALID_WORK_DIMENSION",
		"CL_INVALID_WORK_GROUP_SIZE",
		"CL_INVALID_WORK_ITEM_SIZE",
		"CL_INVALID_GLOBAL_OFFSET",
		"CL_INVALID_EVENT_WAIT_LIST",
		"CL_INVALID_EVENT",
		"CL_INVALID_OPERATION",
		"CL_INVALID_GL_OBJECT",
		"CL_INVALID_BUFFER_SIZE",
		"CL_INVALID_MIP_LEVEL",
		"CL_INVALID_GLOBAL_WORK_SIZE",
		"CL_INVALID_PROPERTY",
		"CL_INVALID_IMAGE_DESCRIPTOR",
		"CL_INVALID_COMPILER_OPTIONS",
		"CL_INVALID_LINKER_OPTIONS",
		"CL_INVALID_DEVICE_PARTITION_COUNT"
	};

	int code = abs(status);
	int err_cnt = sizeof(errors) / sizeof(char*);

	if (code > err_cnt)
		code = 22; //CODE_UNKNOWN

	return errors[code];
}

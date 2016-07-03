#include "astc_codec_batch.h"


#include <stdio.h>

cl_platform_id opencl_platform;
cl_device_id opencl_device;
cl_program opencl_program;
cl_context opencl_context;


static cl_int printDeviceInfo(cl_device_id device)
{
	cl_int status = CL_SUCCESS;
	char buf[512];
	cl_uint clui;
	size_t  clsz;

	status = clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(buf), buf, NULL);
	if (status == CL_SUCCESS) printf("\t\tVendor: %s\n", buf);

	status = clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(buf), buf, NULL);
	if (status == CL_SUCCESS) printf("\t\tName: %s\n", buf);

	status = clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(buf), buf, NULL);
	if (status == CL_SUCCESS) printf("\t\tDriver Version: %s\n", buf);

	status = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &clui, NULL);
	if (status == CL_SUCCESS) printf("\t\tUnits: %i\n", clui);

	status = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &clui, NULL);
	if (status == CL_SUCCESS) printf("\t\tMax dim: %i\n", clui);

	status = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &clsz, NULL);
	if (status == CL_SUCCESS) printf("\t\tWG size: %i\n", clsz);

	status = clGetDeviceInfo(device, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint), &clui, NULL);
	if (status == CL_SUCCESS) printf("\t\tFreqs: %i\n", clui);

	status = clGetDeviceInfo(device, CL_DEVICE_PARTITION_MAX_SUB_DEVICES, sizeof(cl_uint), &clui, NULL);
	if (status == CL_SUCCESS) printf("\t\tMaxPartitions: %i\n", clsz);

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

static cl_int readSource(const char *sourcePath, const char *sourceFilename, char** sourceString)
{
	FILE *fp;
	size_t err;
	size_t size;
#define MAX_PATH_LENGTH 512
	char fullFilename[MAX_PATH_LENGTH];
	strcpy_s(fullFilename, MAX_PATH_LENGTH, sourcePath);
	strcat_s(fullFilename, MAX_PATH_LENGTH, sourceFilename);

	char *source;
	fopen_s(&fp, fullFilename, "rb");
	if (fp == NULL) {
		fprintf(stderr, "Could not open kernel file: %s\n", fullFilename);
		return -1;
	}

	err = fseek(fp, 0, SEEK_END);
	if (err != 0) {
		fprintf(stderr, "Error seeking to end of file\n");
		return -1;
	}

	size = ftell(fp);
	if (size < 0) {
		fprintf(stderr, "Error getting file position\n");
		return -1;
	}

	err = fseek(fp, 0, SEEK_SET);
	if (err != 0) {
		fprintf(stderr, "Error seeking to start of file\n");
		return -1;
	}

	source = new char[size + 1];

	err = fread(source, 1, size, fp);
	if (err != size) {
		fprintf(stderr, "only %d bytes read\n", err);
		delete[] source;
		return -1;
	}

	source[size] = '\0';

	*sourceString = source;
	return CL_SUCCESS;
}

void init_opencl(cl_uint platform_number, cl_uint device_number, int silentmode, int batch_size, int xdim, int ydim, int zdim, int plimit, astc_decode_mode decode_mode)
{
	cl_uint numPlatforms = 0;
	cl_uint numDevices = 0;
	cl_platform_id *platforms;
	cl_device_id *devices;
	cl_int status;

	// get platform
	status = clGetPlatformIDs(0, NULL, &numPlatforms);
	OCL_CHECK_STATUS("Cannot get number of platforms");
	if (numPlatforms < platform_number) {
		fprintf(stderr, "Only %i platforms found. Platform with index %i does not exists.", numPlatforms, platform_number);
		exit(-1);
	}

	platforms = new cl_platform_id[numPlatforms];

	status = clGetPlatformIDs(numPlatforms, platforms, NULL);
	opencl_platform = platforms[platform_number];
	delete[] platforms;

	if (!silentmode) {
		printf("%i platforms found. Using platform %i\n", numPlatforms, platform_number);
		printPlatformInfo(opencl_platform);
	}

	// get device
	status = clGetDeviceIDs(opencl_platform, CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
	OCL_CHECK_STATUS("Cannot get number of devices");
	if (numDevices < device_number) {
		fprintf(stderr, "Only %i devices found. Device with index %i does not exists.", numDevices, device_number);
		exit(-1);
	}

	devices = new cl_device_id[numDevices];

	status = clGetDeviceIDs(opencl_platform, CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);
	opencl_device = devices[device_number];
	delete[] devices;
	if (!silentmode) {
		printf("\n%i devices found on platform %i. Using device %i\n", numDevices, platform_number, device_number);
		printDeviceInfo(opencl_device);
	}

	// create context
	cl_context_properties context_props[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)opencl_platform, 0};
	opencl_context = clCreateContext(context_props, 1, &opencl_device, NULL, NULL, &status);
	OCL_CHECK_STATUS("Cannot create context");

	// create program and compile kernels
	char *source_names[] = { OPENCL_KERNEL_FBP_FILENAME	};
#define FILE_COUNT ((int)(sizeof(source_names)/sizeof(source_names[0])))
	char *sources[FILE_COUNT];
	for (int i = 0; i < FILE_COUNT; i++)
	{
		status = readSource(OPENCL_KERNELS_SOURCE_PATH, source_names[i], &sources[i]);
		OCL_CHECK_STATUS("Cannot read OpenCL kernel source file");
	}

	opencl_program = clCreateProgramWithSource(opencl_context, FILE_COUNT, (const char **)sources, NULL, &status);
	for (int i = 0; i < FILE_COUNT; i++)
		delete[] sources[i];
	OCL_CHECK_STATUS("Cannot create OpenCL program object");


	//set compile-time constants
	int texels_per_block = xdim * ydim * zdim;

	char compile_flags[1024] = "";
	if (sizeof(void*) == 4)
		strcat_s(compile_flags, 1024, " -D OCL_USE_32BIT_POINTERS");

	if (sizeof(void*) == 8)
		strcat_s(compile_flags, 1024, " -D OCL_USE_64BIT_POINTERS");

	char compileOptions[4096];
	sprintf_s(compileOptions, 4096, "%s -I %s -D XDIM=%i -D YDIM=%i -D ZDIM=%i -D TEXELS_PER_BLOCK=%i -D PLIMIT=%i -D FBP_PARTITION_CANDIDATES=%i %s",
		OPENCL_COMPILER_OPTIONS, OPENCL_KERNELS_SOURCE_PATH, xdim, ydim, zdim, texels_per_block, plimit, PARTITION_CANDIDATES, compile_flags);
	if (!silentmode)
		printf("Batch size: %i\nOpenCL compiler options:\n%s\n\n", batch_size, compileOptions);

	//build program
	if (!silentmode)
		printf("Compiling kernels...  ");
	status = clBuildProgram(opencl_program, 1, &opencl_device, compileOptions, NULL, NULL);
	if (status != CL_SUCCESS)
	{
		fprintf(stderr, "Cannot build OpenCL program");

		if (!silentmode) {
			//getting log size
			size_t logSize = 1;
			status = clGetProgramBuildInfo(opencl_program, opencl_device, CL_PROGRAM_BUILD_LOG, NULL, NULL, &logSize);
			OCL_CHECK_STATUS("Cannot get build log size");

			//retreiving build log
			char *buildLog;
			buildLog = new char[logSize];
			
			status = clGetProgramBuildInfo(opencl_program, opencl_device, CL_PROGRAM_BUILD_LOG, logSize, buildLog, NULL);
			if (status != CL_SUCCESS) {
				fprintf(stderr, "Unable to retreive Build log");
				delete[] buildLog;
				exit(-1);
			}
			buildLog[logSize - 1] = '\0';

			fprintf(stderr, "Build Log:\n%s\n", buildLog);
		}

		exit(-1);
	}

	if (!silentmode)
		printf("Done\n\n");
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

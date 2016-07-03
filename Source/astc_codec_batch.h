#ifndef ASTC_CODEC_BATCH_INCLUDED

#define ASTC_CODEC_BATCH_INCLUDED

#include <CL/opencl.h>
#include "astc_codec_internals.h"

#define OPENCL_KERNELS_SOURCE_PATH ".\\"
#define OPENCL_KERNEL_FBP_FILENAME "find_best_partitioning.cl"
#define OPENCL_COMPILER_OPTIONS "-cl-mad-enable"

#define OCL_CHECK_STATUS(str) if(status != CL_SUCCESS) { fprintf(stderr, "%s, errorcode: %i %s\n", str, status, cl_errcode_to_str(status)); exit(-1); }
#define OCL_RELEASE_OBJECT(type, name) { status = clRelease##type(name); OCL_CHECK_STATUS("Cannot release "#type" "#name); }

void init_opencl(cl_uint platform_number, cl_uint device_number, int silentmode, int batch_size, int xdim, int ydim, int zdim, int plimit, astc_decode_mode decode_mode);
void destroy_opencl();
char const* cl_errcode_to_str(cl_int status);

// how many scb candidates will be tested in each compression mode
#define SCB_CANDIDATES 4

// how many partitions will be tested for each multipartition mode
#define PARTITION_CANDIDATES 2

// flags used to skip some compression modes
#define BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF 1
#define BLOCK_STAT_LOWEST_CORREL_CUTOFF (1 << 1)
#define BLOCK_STAT_NORMAL_MAP (1 << 2)
#define BLOCK_STAT_NO_ALPHA (1 << 3)
#define BLOCK_STAT_GREYSCALE (1 << 4)
#define BLOCK_STAT_SKIP_ALL (1 << 5)

// buffers used to store intermediate data in compress_symbolic_block_fixed_partition_*()
struct compress_fixed_partition_buffers
{
	endpoints_and_weights* ei1;
	endpoints_and_weights* ei2;
	endpoints_and_weights* eix1;
	endpoints_and_weights* eix2;
	float *decimated_quantized_weights;
	float *decimated_weights;
	float *flt_quantized_decimated_quantized_weights;
	uint8_t *u8_quantized_decimated_quantized_weights;
};

// buffers and constants used to store intermediate data in find_best_partitionings_batch()
struct find_best_partitionings_buffers
{
	uint16_t * partition_sequence;
	float * uncorr_errors; // partitioning errors assuming uncorrellated-chrominance endpoints
	float * samechroma_errors; // partitioning errors assuming same-chrominance endpoints
	float * separate_errors; // partitioning errors assuming that one of the color channels is uncorrellated from all the other ones

	float weight_imprecision_estim_squared;
	uint16_t partition_search_limits[5];
};


class SymbolicBatchCompressor
{
public:
	SymbolicBatchCompressor(int batch_size, int xdim, int ydim, int zdim, astc_decode_mode decode_mode, const error_weighting_params *ewp);
	void compress_symbolic_batch(const astc_codec_image * input_image, const imageblock * blk_batch, symbolic_compressed_block * scb_batch, int batch_size);
	~SymbolicBatchCompressor();

private:
	const int max_batch_size;
	int batch_size;
	int xdim, ydim, zdim;
	astc_decode_mode decode_mode;
	error_weighting_params ewp;

	uint8_t * blk_stat; //used to skip some compression modes
	
	// buffers to store intermediate data during encoding
	// buffers used in compress_symbolic_batch()
	float * error_weight_sum_batch;
	float * best_errorvals_in_1pl_1partition_mode;
	float * best_errorvals_in_1pl_2partition_mode;
	float * error_of_best_block;
	error_weight_block * ewb_batch;
	symbolic_compressed_block * scb_candidates;
	uint16_t * partition_indices_1plane_batch;
	uint16_t * partition_indices_2planes_batch;

	// buffers used in compress_symbolic_batch_fixed_partition_*()
	compress_fixed_partition_buffers tmpplanes;
	find_best_partitionings_buffers fbp;

	void allocate_buffers(int max_blocks);
	void compress_symbolic_batch_fixed_partition_1_plane(float mode_cutoff, int partition_count, int partition_offset, const imageblock * blk_batch, symbolic_compressed_block * scb_candidates);
	void compress_symbolic_batch_fixed_partition_2_planes(float mode_cutoff, int partition_count, int partition_offset, int separate_component, const imageblock * blk_batch, symbolic_compressed_block * scb_candidates, uint8_t skip_mode);
	void find_best_partitionings_batch(int partition_count, const imageblock * blk_batch, uint16_t * partition_indices_1plane_batch, uint16_t * partition_indices_2planes_batch);
	void find_best_partitionings(int partition_search_limit, int partition_count, const imageblock * pb, const error_weight_block * ewb, uint16_t *best_partitions_single_weight_plane, uint16_t *best_partitions_dual_weight_planes);
};


#endif
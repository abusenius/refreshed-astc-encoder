#ifndef ASTC_CODEC_BATCH_INCLUDED

#define ASTC_CODEC_BATCH_INCLUDED

#include "astc_codec_internals.h"
#include "astc_opencl.h"

// buffers used to store intermediate data in compress_symbolic_block_fixed_partition_*()
struct compress_fixed_partition_buffers
{
	endpoints_and_weights* ei1;
	endpoints_and_weights* ei2;
	endpoints_and_weights* eix1;
	endpoints_and_weights* eix2;
	float *flt_quantized_decimated_quantized_weights;
	uint8_t *u8_quantized_decimated_quantized_weights;
	uint8_t *scb_stat;
	struct
	{
		endpoints* ep;
		int *partition_format_specifiers;
		int *weight_mode;
		int *color_quantization_level;
		int *color_quantization_level_mod;
		int *decimation_mode;
		uint8_t *u8_qdq_weights;
	} per_scb;
};

// buffers and constants used to store intermediate data in find_best_partitionings_batch()
struct find_best_partitionings_buffers
{
	cl_kernel find_best_partitionings;

	ocl_buffer<uint16_t, ocl_buffer_type::DEVICE> partition_sequence;
	cl_mem ptab[5];
	cl_mem ptab_pseudo[3];

	float weight_imprecision_estim_squared;
	uint16_t partition_search_limits[5];
};

struct compute_angular_endpoints_buffers
{
	cl_kernel compute_angular_endpoints_1plane;
	cl_kernel compute_angular_endpoints_2planes;
	cl_mem bsd_1plane;
	cl_mem bsd_2planes;
	ocl_buffer<float, ocl_buffer_type::DEVICE> decimated_quantized_weights;
	ocl_buffer<float, ocl_buffer_type::DEVICE> decimated_weights;
	ocl_buffer<float, ocl_buffer_type::DEVICE> weight_low_value1;
	ocl_buffer<float, ocl_buffer_type::DEVICE> weight_low_value2;
	ocl_buffer<float, ocl_buffer_type::DEVICE> weight_high_value1;
	ocl_buffer<float, ocl_buffer_type::DEVICE> weight_high_value2;
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

	cl_command_queue opencl_queue;

	ocl_buffer<uint8_t, ocl_buffer_type::DEVICE> blk_stat; //used to skip some compression modes
	cl_mem blk_buf;
	
	// buffers to store intermediate data during encoding
	// buffers used in compress_symbolic_batch()
	float * error_weight_sum_batch;
	float * best_errorvals_in_1pl_1partition_mode;
	float * best_errorvals_in_1pl_2partition_mode;
	float * error_of_best_block;
	ocl_buffer<error_weight_block, ocl_buffer_type::DEVICE> ewb_batch;
	symbolic_compressed_block * scb_candidates;
	ocl_buffer<uint16_t, ocl_buffer_type::DEVICE> partition_indices_1plane_batch;
	ocl_buffer<uint16_t, ocl_buffer_type::DEVICE> partition_indices_2planes_batch;

	//ocl_buffer<int4, ocl_buffer_type::DEVICE> idebug;
	//ocl_buffer<float4, ocl_buffer_type::DEVICE> fdebug;
	void compute_angular_endpoints_1plane_batch_ocl();
	void compute_angular_endpoints_2planes_batch_ocl(uint8_t skip_mode);

	// buffers used in compress_symbolic_batch_fixed_partition_*()
	compress_fixed_partition_buffers tmpplanes;
	find_best_partitionings_buffers fbp;
	compute_angular_endpoints_buffers cae;

	void allocate_buffers(int max_blocks);
	void compress_symbolic_batch_fixed_partition_1_plane(int partition_count, int empty_partition_count, int partition_offset, const imageblock * blk_batch, symbolic_compressed_block * scb_candidates);
	void compress_symbolic_batch_fixed_partition_2_planes(int partition_count, int empty_partition_count, int partition_offset, int separate_component, const imageblock * blk_batch, symbolic_compressed_block * scb_candidates, uint8_t skip_mode);
	void find_best_partitionings_batch(int partition_count, int empty_partition_count, const imageblock * blk_batch);
	void find_best_partitionings_batch_ocl(int partition_count, int empty_partition_count, const imageblock * blk_batch);
	void find_best_partitionings(int partition_search_limit, int empty_partition_count, int partition_count, const imageblock * pb, const error_weight_block * ewb, uint16_t *best_partitions_single_weight_plane, uint16_t *best_partitions_dual_weight_planes);
};


#endif
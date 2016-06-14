#ifndef ASTC_CODEC_BATCH_INCLUDED

#define ASTC_CODEC_BATCH_INCLUDED


#include "astc_codec_internals.h"

// how many scb candidates will be tested in each compression mode
#define SCB_CANDIDATES 4

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



class SymbolicBatchCompressor
{
public:
	SymbolicBatchCompressor(int batch_size);
	SymbolicBatchCompressor(int batch_size, int xdim, int ydim, int zdim, astc_decode_mode decode_mode, const error_weighting_params *ewp);
	//void set_tile_size(int x, int y, int z) { xdim = x; ydim = y; zdim = z; };
	//void set_decode_mode(astc_decode_mode mode) { decode_mode = mode; };
	//void set_error_weighting_params(const error_weighting_params *params) { ewp = *params; };
	void compress_symbolic_batch(const astc_codec_image * input_image, const imageblock * blk_batch, symbolic_compressed_block * scb_batch, int batch_size);
	~SymbolicBatchCompressor();

private:
	const int max_batch_size;
	int xdim, ydim, zdim;
	astc_decode_mode decode_mode;
	error_weighting_params ewp;

	uint8_t * blk_finished; // 1 if block is compressed well enough
	uint8_t * blk_skip_2planes; // 1 if decided not to test 2 planes mode or *blk_finished==1

	// buffers to store intermediate data during encoding
	// buffers used in compress_symbolic_batch()
	error_weight_block * ewb_batch;
	symbolic_compressed_block * scb_candidates;

	// buffers used in compress_symbolic_batch_fixed_partition_*()
	compress_fixed_partition_buffers tmpplanes;

	void allocate_buffers(int max_blocks);
};


#endif
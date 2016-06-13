#ifndef ASTC_CODEC_BATCH_INCLUDED

#define ASTC_CODEC_BATCH_INCLUDED


#include "astc_codec_internals.h"

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

struct compress_symbolic_block_buffers
{
	error_weight_block *ewb;
	symbolic_compressed_block *tempblocks;
	imageblock *temp;
	compress_fixed_partition_buffers *plane1;
	compress_fixed_partition_buffers *planes2;
};


float compress_symbolic_block(const astc_codec_image * input_image,
	astc_decode_mode decode_mode, int xdim, int ydim, int zdim, const error_weighting_params * ewp, const imageblock * blk, symbolic_compressed_block * scb,
	compress_symbolic_block_buffers * tmpbuf);

class SymbolicBatchCompressor
{
public:
	SymbolicBatchCompressor(int batch_size);
	SymbolicBatchCompressor(int batch_size, int xdim, int ydim, int zdim, astc_decode_mode decode_mode, const error_weighting_params *ewp);
	//void set_tile_size(int x, int y, int z) { xdim = x; ydim = y; zdim = z; };
	//void set_decode_mode(astc_decode_mode mode) { decode_mode = mode; };
	//void set_error_weighting_params(const error_weighting_params *params) { ewp = *params; };
	void compress_symbolic_batch(const astc_codec_image * input_image, const imageblock * blk, symbolic_compressed_block * scb, int batch_size);
	~SymbolicBatchCompressor();

private:
	const int max_batch_size;
	int xdim, ydim, zdim;
	astc_decode_mode decode_mode;
	error_weighting_params ewp;

	//buffers to store intermediate data during encoding
	compress_symbolic_block_buffers tmpbuf;
	compress_fixed_partition_buffers tmpplanes;

	void allocate_buffers(int max_blocks);
};


#endif
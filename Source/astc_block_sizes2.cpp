/*----------------------------------------------------------------------------*/  
/**
 *	This confidential and proprietary software may be used only as
 *	authorised by a licensing agreement from ARM Limited
 *	(C) COPYRIGHT 2011-2012 ARM Limited
 *	ALL RIGHTS RESERVED
 *
 *	The entire notice above must be reproduced on all authorised
 *	copies and copies may only be made to the extent permitted
 *	by a licensing agreement from ARM Limited.
 *
 *	@brief	For ASTC, generate the block size descriptor and the associated
 *			decimation tables.
 */ 
/*----------------------------------------------------------------------------*/ 

#include <assert.h>
#include "astc_codec_internals.h"
#include "math.h"

struct block_size_info
{
	block_size_descriptor bsd;
	block_size_descriptor_sorted bsd_1plane;
	block_size_descriptor_sorted bsd_2planes;
};

static block_size_info *bsd_pointers[4096];

extern const float percentile_table_4x4[2048];
extern const float percentile_table_4x5[2048];
extern const float percentile_table_4x6[2048];
extern const float percentile_table_4x8[2048];
extern const float percentile_table_4x10[2048];
extern const float percentile_table_4x12[2048];
extern const float percentile_table_5x4[2048];
extern const float percentile_table_5x5[2048];
extern const float percentile_table_5x6[2048];
extern const float percentile_table_5x8[2048];
extern const float percentile_table_5x10[2048];
extern const float percentile_table_5x12[2048];
extern const float percentile_table_6x4[2048];
extern const float percentile_table_6x5[2048];
extern const float percentile_table_6x6[2048];
extern const float percentile_table_6x8[2048];
extern const float percentile_table_6x10[2048];
extern const float percentile_table_6x12[2048];
extern const float percentile_table_8x4[2048];
extern const float percentile_table_8x5[2048];
extern const float percentile_table_8x6[2048];
extern const float percentile_table_8x8[2048];
extern const float percentile_table_8x10[2048];
extern const float percentile_table_8x12[2048];
extern const float percentile_table_10x4[2048];
extern const float percentile_table_10x5[2048];
extern const float percentile_table_10x6[2048];
extern const float percentile_table_10x8[2048];
extern const float percentile_table_10x10[2048];
extern const float percentile_table_10x12[2048];
extern const float percentile_table_12x4[2048];
extern const float percentile_table_12x5[2048];
extern const float percentile_table_12x6[2048];
extern const float percentile_table_12x8[2048];
extern const float percentile_table_12x10[2048];
extern const float percentile_table_12x12[2048];

const float *get_2d_percentile_table(int blockdim_x, int blockdim_y)
{
	switch (blockdim_x)
	{
	case 4:
		switch (blockdim_y)
		{
		case 4:
			return percentile_table_4x4;
		case 5:
			return percentile_table_4x5;
		case 6:
			return percentile_table_4x6;
		case 8:
			return percentile_table_4x8;
		case 10:
			return percentile_table_4x10;
		case 12:
			return percentile_table_4x12;
		}
		break;
	case 5:
		switch (blockdim_y)
		{
		case 4:
			return percentile_table_5x4;
		case 5:
			return percentile_table_5x5;
		case 6:
			return percentile_table_5x6;
		case 8:
			return percentile_table_5x8;
		case 10:
			return percentile_table_5x10;
		case 12:
			return percentile_table_5x12;
		}
		break;

	case 6:
		switch (blockdim_y)
		{
		case 4:
			return percentile_table_6x4;
		case 5:
			return percentile_table_6x5;
		case 6:
			return percentile_table_6x6;
		case 8:
			return percentile_table_6x8;
		case 10:
			return percentile_table_6x10;
		case 12:
			return percentile_table_6x12;
		}
		break;

	case 8:
		switch (blockdim_y)
		{
		case 4:
			return percentile_table_8x4;
		case 5:
			return percentile_table_8x5;
		case 6:
			return percentile_table_8x6;
		case 8:
			return percentile_table_8x8;
		case 10:
			return percentile_table_8x10;
		case 12:
			return percentile_table_8x12;
		}
		break;

	case 10:
		switch (blockdim_y)
		{
		case 4:
			return percentile_table_10x4;
		case 5:
			return percentile_table_10x5;
		case 6:
			return percentile_table_10x6;
		case 8:
			return percentile_table_10x8;
		case 10:
			return percentile_table_10x10;
		case 12:
			return percentile_table_10x12;
		}
		break;

	case 12:
		switch (blockdim_y)
		{
		case 4:
			return percentile_table_12x4;
		case 5:
			return percentile_table_12x5;
		case 6:
			return percentile_table_12x6;
		case 8:
			return percentile_table_12x8;
		case 10:
			return percentile_table_12x10;
		case 12:
			return percentile_table_12x12;
		}
		break;
	default:
		break;
	}

	return NULL;				// should never happen.
}

// stubbed for the time being.
static const float dummy_percentile_table_3d[2048] = { 0 };
const float *get_3d_percentile_table(int blockdim_x, int blockdim_y, int blockdim_z)
{
	IGNORE(blockdim_x);
	IGNORE(blockdim_y);
	IGNORE(blockdim_z);
	return dummy_percentile_table_3d;
}



// return 0 on invalid mode, 1 on valid mode.
static int decode_block_mode_2d(int blockmode, int *Nval, int *Mval, int *dual_weight_plane, int *quant_mode)
{
	int base_quant_mode = (blockmode >> 4) & 1;
	int H = (blockmode >> 9) & 1;
	int D = (blockmode >> 10) & 1;

	int A = (blockmode >> 5) & 0x3;

	int N = 0, M = 0;

	if ((blockmode & 3) != 0)
	{
		base_quant_mode |= (blockmode & 3) << 1;
		int B = (blockmode >> 7) & 3;
		switch ((blockmode >> 2) & 3)
		{
		case 0:
			N = B + 4;
			M = A + 2;
			break;
		case 1:
			N = B + 8;
			M = A + 2;
			break;
		case 2:
			N = A + 2;
			M = B + 8;
			break;
		case 3:
			B &= 1;
			if (blockmode & 0x100)
			{
				N = B + 2;
				M = A + 2;
			}
			else
			{
				N = A + 2;
				M = B + 6;
			}
			break;
		}
	}
	else
	{
		base_quant_mode |= ((blockmode >> 2) & 3) << 1;
		if (((blockmode >> 2) & 3) == 0)
			return 0;
		int B = (blockmode >> 9) & 3;
		switch ((blockmode >> 7) & 3)
		{
		case 0:
			N = 12;
			M = A + 2;
			break;
		case 1:
			N = A + 2;
			M = 12;
			break;
		case 2:
			N = A + 6;
			M = B + 6;
			D = 0;
			H = 0;
			break;
		case 3:
			switch ((blockmode >> 5) & 3)
			{
			case 0:
				N = 6;
				M = 10;
				break;
			case 1:
				N = 10;
				M = 6;
				break;
			case 2:
			case 3:
				return 0;
			}
			break;
		}
	}

	int weight_count = N * M * (D + 1);
	int qmode = (base_quant_mode - 2) + 6 * H;

	int weightbits = compute_ise_bitcount(weight_count, (quantization_method) qmode);
	if (weight_count > MAX_WEIGHTS_PER_BLOCK || weightbits < MIN_WEIGHT_BITS_PER_BLOCK || weightbits > MAX_WEIGHT_BITS_PER_BLOCK)
		return 0;

	*Nval = N;
	*Mval = M;
	*dual_weight_plane = D;
	*quant_mode = qmode;
	return 1;
}


static int decode_block_mode_3d(int blockmode, int *Nval, int *Mval, int *Qval, int *dual_weight_plane, int *quant_mode)
{
	int base_quant_mode = (blockmode >> 4) & 1;
	int H = (blockmode >> 9) & 1;
	int D = (blockmode >> 10) & 1;

	int A = (blockmode >> 5) & 0x3;

	int N = 0, M = 0, Q = 0;

	if ((blockmode & 3) != 0)
	{
		base_quant_mode |= (blockmode & 3) << 1;
		int B = (blockmode >> 7) & 3;
		int C = (blockmode >> 2) & 0x3;
		N = A + 2;
		M = B + 2;
		Q = C + 2;
	}
	else
	{
		base_quant_mode |= ((blockmode >> 2) & 3) << 1;
		if (((blockmode >> 2) & 3) == 0)
			return 0;
		int B = (blockmode >> 9) & 3;
		if (((blockmode >> 7) & 3) != 3)
		{
			D = 0;
			H = 0;
		}
		switch ((blockmode >> 7) & 3)
		{
		case 0:
			N = 6;
			M = B + 2;
			Q = A + 2;
			break;
		case 1:
			N = A + 2;
			M = 6;
			Q = B + 2;
			break;
		case 2:
			N = A + 2;
			M = B + 2;
			Q = 6;
			break;
		case 3:
			N = 2;
			M = 2;
			Q = 2;
			switch ((blockmode >> 5) & 3)
			{
			case 0:
				N = 6;
				break;
			case 1:
				M = 6;
				break;
			case 2:
				Q = 6;
				break;
			case 3:
				return 0;
			}
			break;
		}
	}

	int weight_count = N * M * Q * (D + 1);
	int qmode = (base_quant_mode - 2) + 6 * H;

	int weightbits = compute_ise_bitcount(weight_count, (quantization_method) qmode);
	if (weight_count > MAX_WEIGHTS_PER_BLOCK || weightbits < MIN_WEIGHT_BITS_PER_BLOCK || weightbits > MAX_WEIGHT_BITS_PER_BLOCK)
		return 0;

	*Nval = N;
	*Mval = M;
	*Qval = Q;
	*dual_weight_plane = D;
	*quant_mode = qmode;
	return 1;
}




static void initialize_decimation_table_2d(
											  // dimensions of the block
											  int xdim, int ydim,
											  // number of grid points in 2d weight grid
											  int x_weights, int y_weights, decimation_table * dt)
{
	int i, j;
	int x, y;

	int texels_per_block = xdim * ydim;
	int weights_per_block = x_weights * y_weights;

	int weightcount_of_texel[MAX_TEXELS_PER_BLOCK];
	int grid_weights_of_texel[MAX_TEXELS_PER_BLOCK][4];
	int weights_of_texel[MAX_TEXELS_PER_BLOCK][4];

	int texelcount_of_weight[MAX_WEIGHTS_PER_BLOCK];
	int texels_of_weight[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK];
	int texelweights_of_weight[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK];

	for (i = 0; i < weights_per_block; i++)
		texelcount_of_weight[i] = 0;
	for (i = 0; i < texels_per_block; i++)
		weightcount_of_texel[i] = 0;

	for (y = 0; y < ydim; y++)
		for (x = 0; x < xdim; x++)
		{
			int texel = y * xdim + x;

			int x_weight = (((1024 + xdim / 2) / (xdim - 1)) * x * (x_weights - 1) + 32) >> 6;
			int y_weight = (((1024 + ydim / 2) / (ydim - 1)) * y * (y_weights - 1) + 32) >> 6;

			int x_weight_frac = x_weight & 0xF;
			int y_weight_frac = y_weight & 0xF;
			int x_weight_int = x_weight >> 4;
			int y_weight_int = y_weight >> 4;
			int qweight[4];
			int weight[4];
			qweight[0] = x_weight_int + y_weight_int * x_weights;
			qweight[1] = qweight[0] + 1;
			qweight[2] = qweight[0] + x_weights;
			qweight[3] = qweight[2] + 1;

			// truncated-precision bilinear interpolation.
			int prod = x_weight_frac * y_weight_frac;

			weight[3] = (prod + 8) >> 4;
			weight[1] = x_weight_frac - weight[3];
			weight[2] = y_weight_frac - weight[3];
			weight[0] = 16 - x_weight_frac - y_weight_frac + weight[3];

			for (i = 0; i < 4; i++)
				if (weight[i] != 0)
				{
					grid_weights_of_texel[texel][weightcount_of_texel[texel]] = qweight[i];
					weights_of_texel[texel][weightcount_of_texel[texel]] = weight[i];
					weightcount_of_texel[texel]++;
					texels_of_weight[qweight[i]][texelcount_of_weight[qweight[i]]] = texel;
					texelweights_of_weight[qweight[i]][texelcount_of_weight[qweight[i]]] = weight[i];
					texelcount_of_weight[qweight[i]]++;
				}
		}

	for (i = 0; i < texels_per_block; i++)
	{
		dt->texel_num_weights[i] = weightcount_of_texel[i];

		// ensure that all 4 entries are actually initialized.
		// This allows a branch-free implemntation of compute_value_of_texel_flt()
		for (j = 0; j < 4; j++)
		{
			dt->texel_weights_int[i][j] = 0;
			dt->texel_weights_float[i][j] = 0.0f;
			dt->texel_weights[i][j] = 0;
		}

		for (j = 0; j < weightcount_of_texel[i]; j++)
		{
			dt->texel_weights_int[i][j] = weights_of_texel[i][j];
			dt->texel_weights_float[i][j] = static_cast < float >(weights_of_texel[i][j]) * (1.0f / TEXEL_WEIGHT_SUM);
			dt->texel_weights[i][j] = grid_weights_of_texel[i][j];
		}
	}

	for (i = 0; i < weights_per_block; i++)
	{
		dt->weight_num_texels[i] = texelcount_of_weight[i];


		for (j = 0; j < texelcount_of_weight[i]; j++)
		{
			dt->weight_texel[i][j] = texels_of_weight[i][j];
			dt->weights_int[i][j] = texelweights_of_weight[i][j];
			dt->weights_flt[i][j] = static_cast < float >(texelweights_of_weight[i][j]);
		}
	}

	dt->num_texels = texels_per_block;
	dt->num_weights = weights_per_block;


}




static void initialize_decimation_table_3d(
											  // dimensions of the block
											  int xdim, int ydim, int zdim,
											  // number of grid points in 3d weight grid
											  int x_weights, int y_weights, int z_weights, decimation_table * dt)
{
	int i, j;
	int x, y, z;

	int texels_per_block = xdim * ydim * zdim;
	int weights_per_block = x_weights * y_weights * z_weights;

	int weightcount_of_texel[MAX_TEXELS_PER_BLOCK];
	int grid_weights_of_texel[MAX_TEXELS_PER_BLOCK][4];
	int weights_of_texel[MAX_TEXELS_PER_BLOCK][4];

	int texelcount_of_weight[MAX_WEIGHTS_PER_BLOCK];
	int texels_of_weight[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK];
	int texelweights_of_weight[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK];

	for (i = 0; i < weights_per_block; i++)
		texelcount_of_weight[i] = 0;
	for (i = 0; i < texels_per_block; i++)
		weightcount_of_texel[i] = 0;

	for (z = 0; z < zdim; z++)
		for (y = 0; y < ydim; y++)
			for (x = 0; x < xdim; x++)
			{
				int texel = (z * ydim + y) * xdim + x;

				int x_weight = (((1024 + xdim / 2) / (xdim - 1)) * x * (x_weights - 1) + 32) >> 6;
				int y_weight = (((1024 + ydim / 2) / (ydim - 1)) * y * (y_weights - 1) + 32) >> 6;
				int z_weight = (((1024 + zdim / 2) / (zdim - 1)) * z * (z_weights - 1) + 32) >> 6;

				int x_weight_frac = x_weight & 0xF;
				int y_weight_frac = y_weight & 0xF;
				int z_weight_frac = z_weight & 0xF;
				int x_weight_int = x_weight >> 4;
				int y_weight_int = y_weight >> 4;
				int z_weight_int = z_weight >> 4;
				int qweight[4];
				int weight[4];
				qweight[0] = (z_weight_int * y_weights + y_weight_int) * x_weights + x_weight_int;
				qweight[3] = ((z_weight_int + 1) * y_weights + (y_weight_int + 1)) * x_weights + (x_weight_int + 1);

				// simplex interpolation
				int fs = x_weight_frac;
				int ft = y_weight_frac;
				int fp = z_weight_frac;

				int cas = ((fs > ft) << 2) + ((ft > fp) << 1) + ((fs > fp));
				int N = x_weights;
				int NM = x_weights * y_weights;

				int s1, s2, w0, w1, w2, w3;
				switch (cas)
				{
				case 7:
					s1 = 1;
					s2 = N;
					w0 = 16 - fs;
					w1 = fs - ft;
					w2 = ft - fp;
					w3 = fp;
					break;
				case 3:
					s1 = N;
					s2 = 1;
					w0 = 16 - ft;
					w1 = ft - fs;
					w2 = fs - fp;
					w3 = fp;
					break;
				case 5:
					s1 = 1;
					s2 = NM;
					w0 = 16 - fs;
					w1 = fs - fp;
					w2 = fp - ft;
					w3 = ft;
					break;
				case 4:
					s1 = NM;
					s2 = 1;
					w0 = 16 - fp;
					w1 = fp - fs;
					w2 = fs - ft;
					w3 = ft;
					break;
				case 2:
					s1 = N;
					s2 = NM;
					w0 = 16 - ft;
					w1 = ft - fp;
					w2 = fp - fs;
					w3 = fs;
					break;
				case 0:
					s1 = NM;
					s2 = N;
					w0 = 16 - fp;
					w1 = fp - ft;
					w2 = ft - fs;
					w3 = fs;
					break;

				default:
					s1 = NM;
					s2 = N;
					w0 = 16 - fp;
					w1 = fp - ft;
					w2 = ft - fs;
					w3 = fs;
					break;
				}

				qweight[1] = qweight[0] + s1;
				qweight[2] = qweight[1] + s2;
				weight[0] = w0;
				weight[1] = w1;
				weight[2] = w2;
				weight[3] = w3;

				/* 
				   for(i=0;i<4;i++) weight[i] <<= 4; */

				for (i = 0; i < 4; i++)
					if (weight[i] != 0)
					{
						grid_weights_of_texel[texel][weightcount_of_texel[texel]] = qweight[i];
						weights_of_texel[texel][weightcount_of_texel[texel]] = weight[i];
						weightcount_of_texel[texel]++;
						texels_of_weight[qweight[i]][texelcount_of_weight[qweight[i]]] = texel;
						texelweights_of_weight[qweight[i]][texelcount_of_weight[qweight[i]]] = weight[i];
						texelcount_of_weight[qweight[i]]++;
					}
			}

	for (i = 0; i < texels_per_block; i++)
	{
		dt->texel_num_weights[i] = weightcount_of_texel[i];

		// ensure that all 4 entries are actually initialized.
		// This allows a branch-free implemntation of compute_value_of_texel_flt()
		for (j = 0; j < 4; j++)
		{
			dt->texel_weights_int[i][j] = 0;
			dt->texel_weights_float[i][j] = 0.0f;
			dt->texel_weights[i][j] = 0;
		}

		for (j = 0; j < weightcount_of_texel[i]; j++)
		{
			dt->texel_weights_int[i][j] = weights_of_texel[i][j];
			dt->texel_weights_float[i][j] = static_cast < float >(weights_of_texel[i][j]) * (1.0f / TEXEL_WEIGHT_SUM);
			dt->texel_weights[i][j] = grid_weights_of_texel[i][j];
		}
	}

	for (i = 0; i < weights_per_block; i++)
	{
		dt->weight_num_texels[i] = texelcount_of_weight[i];
		for (j = 0; j < texelcount_of_weight[i]; j++)
		{
			dt->weight_texel[i][j] = texels_of_weight[i][j];
			dt->weights_int[i][j] = texelweights_of_weight[i][j];
			dt->weights_flt[i][j] = static_cast < float >(texelweights_of_weight[i][j]);
		}
	}

	dt->num_texels = texels_per_block;
	dt->num_weights = weights_per_block;
}



void construct_block_size_descriptor_2d(int xdim, int ydim, block_size_descriptor * bsd)
{
	int decimation_mode_index[256];	// for each of the 256 entries in the decim_table_array, its index
	int decimation_mode_count = 0;

	int i;
	int x_weights;
	int y_weights;

	for (i = 0; i < 256; i++)
	{
		decimation_mode_index[i] = -1;
	}

	// gather all the infill-modes that can be used with the current block size
	for (x_weights = 2; x_weights <= 12; x_weights++)
		for (y_weights = 2; y_weights <= 12; y_weights++)
		{
			if (x_weights * y_weights > MAX_WEIGHTS_PER_BLOCK)
				continue;
			decimation_table *dt = new decimation_table;
			decimation_mode_index[y_weights * 16 + x_weights] = decimation_mode_count;
			initialize_decimation_table_2d(xdim, ydim, x_weights, y_weights, dt);

			int weight_count = x_weights * y_weights;

			int maxprec_1plane = -1;
			int maxprec_2planes = -1;
			for (i = 0; i < 12; i++)
			{
				int bits_1plane = compute_ise_bitcount(weight_count, (quantization_method) i);
				int bits_2planes = compute_ise_bitcount(2 * weight_count, (quantization_method) i);
				if (bits_1plane >= MIN_WEIGHT_BITS_PER_BLOCK && bits_1plane <= MAX_WEIGHT_BITS_PER_BLOCK)
					maxprec_1plane = i;
				if (bits_2planes >= MIN_WEIGHT_BITS_PER_BLOCK && bits_2planes <= MAX_WEIGHT_BITS_PER_BLOCK)
					maxprec_2planes = i;
			}

			if (2 * x_weights * y_weights > MAX_WEIGHTS_PER_BLOCK)
				maxprec_2planes = -1;

			bsd->permit_encode[decimation_mode_count] = (x_weights <= xdim && y_weights <= ydim);

			bsd->decimation_mode_samples[decimation_mode_count] = weight_count;
			bsd->decimation_mode_maxprec_1plane[decimation_mode_count] = maxprec_1plane;
			bsd->decimation_mode_maxprec_2planes[decimation_mode_count] = maxprec_2planes;
			bsd->decimation_tables[decimation_mode_count] = dt;

			decimation_mode_count++;
		}

	for (i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		bsd->decimation_mode_percentile_1plane[i] = INFINITY;
		bsd->decimation_mode_percentile_2planes[i] = INFINITY;
	}

	for (i = decimation_mode_count; i < MAX_DECIMATION_MODES; i++)
	{
		bsd->permit_encode[i] = 0;
		bsd->decimation_mode_samples[i] = 0;
		bsd->decimation_mode_maxprec_1plane[i] = -1;
		bsd->decimation_mode_maxprec_2planes[i] = -1;
	}

	bsd->decimation_mode_count = decimation_mode_count;

	const float *percentiles = get_2d_percentile_table(xdim, ydim);

	// then construct the list of block formats
	for (i = 0; i < 2048; i++)
	{
		int x_weights, y_weights;
		int is_dual_plane;
		int quantization_mode;
		int fail = 0;
		int permit_encode = 1;

		if (decode_block_mode_2d(i, &x_weights, &y_weights, &is_dual_plane, &quantization_mode))
		{
			if (x_weights > xdim || y_weights > ydim)
				permit_encode = 0;
		}
		else
		{
			fail = 1;
			permit_encode = 0;
		}

		if (fail)
		{
			bsd->block_modes[i].decimation_mode = -1;
			bsd->block_modes[i].quantization_mode = -1;
			bsd->block_modes[i].is_dual_plane = -1;
			bsd->block_modes[i].permit_encode = 0;
			bsd->block_modes[i].permit_decode = 0;
			bsd->block_modes[i].percentile = 1.0f;
		}
		else
		{
			int decimation_mode = decimation_mode_index[y_weights * 16 + x_weights];
			bsd->block_modes[i].decimation_mode = decimation_mode;
			bsd->block_modes[i].quantization_mode = quantization_mode;
			bsd->block_modes[i].is_dual_plane = is_dual_plane;
			bsd->block_modes[i].permit_encode = permit_encode;
			bsd->block_modes[i].permit_decode = permit_encode;	// disallow decode of grid size larger than block size.
			bsd->block_modes[i].percentile = percentiles[i];

			if (bsd->decimation_mode_percentile_1plane[decimation_mode] > percentiles[i] && !is_dual_plane && permit_encode)
				bsd->decimation_mode_percentile_1plane[decimation_mode] = percentiles[i];

			if (bsd->decimation_mode_percentile_2planes[decimation_mode] > percentiles[i] && is_dual_plane && permit_encode)
				bsd->decimation_mode_percentile_2planes[decimation_mode] = percentiles[i];
		}

	}

	if (xdim * ydim <= 64)
	{
		bsd->texelcount_for_bitmap_partitioning = xdim * ydim;
		for (i = 0; i < xdim * ydim; i++)
			bsd->texels_for_bitmap_partitioning[i] = i;
	}

	else
	{
		// pick 64 random texels for use with bitmap partitioning.
		int arr[MAX_TEXELS_PER_BLOCK];
		for (i = 0; i < xdim * ydim; i++)
			arr[i] = 0;
		int arr_elements_set = 0;
		while (arr_elements_set < 64)
		{
			int idx = rand() % (xdim * ydim);
			if (arr[idx] == 0)
			{
				arr_elements_set++;
				arr[idx] = 1;
			}
		}
		int texel_weights_written = 0;
		int idx = 0;
		while (texel_weights_written < 64)
		{
			if (arr[idx])
				bsd->texels_for_bitmap_partitioning[texel_weights_written++] = idx;
			idx++;
		}
		bsd->texelcount_for_bitmap_partitioning = 64;

	}
}



void construct_block_size_descriptor_3d(int xdim, int ydim, int zdim, block_size_descriptor * bsd)
{
	int decimation_mode_index[512];	// for each of the 512 entries in the decim_table_array, its index
	int decimation_mode_count = 0;

	int i;
	int x_weights;
	int y_weights;
	int z_weights;

	for (i = 0; i < 512; i++)
	{
		decimation_mode_index[i] = -1;
	}

	// gather all the infill-modes that can be used with the current block size
	for (x_weights = 2; x_weights <= 6; x_weights++)
		for (y_weights = 2; y_weights <= 6; y_weights++)
			for (z_weights = 2; z_weights <= 6; z_weights++)
			{
				if ((x_weights * y_weights * z_weights) > MAX_WEIGHTS_PER_BLOCK)
					continue;
				decimation_table *dt = new decimation_table;
				decimation_mode_index[z_weights * 64 + y_weights * 8 + x_weights] = decimation_mode_count;
				initialize_decimation_table_3d(xdim, ydim, zdim, x_weights, y_weights, z_weights, dt);

				int weight_count = x_weights * y_weights * z_weights;

				int maxprec_1plane = -1;
				int maxprec_2planes = -1;
				for (i = 0; i < 12; i++)
				{
					int bits_1plane = compute_ise_bitcount(weight_count, (quantization_method) i);
					int bits_2planes = compute_ise_bitcount(2 * weight_count, (quantization_method) i);
					if (bits_1plane >= MIN_WEIGHT_BITS_PER_BLOCK && bits_1plane <= MAX_WEIGHT_BITS_PER_BLOCK)
						maxprec_1plane = i;
					if (bits_2planes >= MIN_WEIGHT_BITS_PER_BLOCK && bits_2planes <= MAX_WEIGHT_BITS_PER_BLOCK)
						maxprec_2planes = i;
				}

				if ((2 * x_weights * y_weights * z_weights) > MAX_WEIGHTS_PER_BLOCK)
					maxprec_2planes = -1;

				bsd->permit_encode[decimation_mode_count] = (x_weights <= xdim && y_weights <= ydim && z_weights <= zdim);

				bsd->decimation_mode_samples[decimation_mode_count] = weight_count;
				bsd->decimation_mode_maxprec_1plane[decimation_mode_count] = maxprec_1plane;
				bsd->decimation_mode_maxprec_2planes[decimation_mode_count] = maxprec_2planes;
				bsd->decimation_tables[decimation_mode_count] = dt;

				decimation_mode_count++;
			}

	for (i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		bsd->decimation_mode_percentile_1plane[i] = INFINITY;
		bsd->decimation_mode_percentile_2planes[i] = INFINITY;
	}

	for (i = decimation_mode_count; i < MAX_DECIMATION_MODES; i++)
	{
		bsd->permit_encode[i] = 0;
		bsd->decimation_mode_samples[i] = 0;
		bsd->decimation_mode_maxprec_1plane[i] = -1;
		bsd->decimation_mode_maxprec_2planes[i] = -1;
	}

	bsd->decimation_mode_count = decimation_mode_count;

	const float *percentiles = get_3d_percentile_table(xdim, ydim, zdim);

	// then construct the list of block formats
	for (i = 0; i < 2048; i++)
	{
		int x_weights, y_weights, z_weights;
		int is_dual_plane;
		int quantization_mode;
		int fail = 0;
		int permit_encode = 1;

		if (decode_block_mode_3d(i, &x_weights, &y_weights, &z_weights, &is_dual_plane, &quantization_mode))
		{
			if (x_weights > xdim || y_weights > ydim || z_weights > zdim)
				permit_encode = 0;
		}
		else
		{
			fail = 1;
			permit_encode = 0;
		}
		if (fail)
		{
			bsd->block_modes[i].decimation_mode = -1;
			bsd->block_modes[i].quantization_mode = -1;
			bsd->block_modes[i].is_dual_plane = -1;
			bsd->block_modes[i].permit_encode = 0;
			bsd->block_modes[i].permit_decode = 0;
			bsd->block_modes[i].percentile = 1.0f;
		}
		else
		{
			int decimation_mode = decimation_mode_index[z_weights * 64 + y_weights * 8 + x_weights];
			bsd->block_modes[i].decimation_mode = decimation_mode;
			bsd->block_modes[i].quantization_mode = quantization_mode;
			bsd->block_modes[i].is_dual_plane = is_dual_plane;
			bsd->block_modes[i].permit_encode = permit_encode;
			bsd->block_modes[i].permit_decode = permit_encode;
			bsd->block_modes[i].percentile = percentiles[i];

			if (bsd->decimation_mode_percentile_1plane[decimation_mode] > percentiles[i] && !is_dual_plane && permit_encode)
				bsd->decimation_mode_percentile_1plane[decimation_mode] = percentiles[i];

			if (bsd->decimation_mode_percentile_2planes[decimation_mode] > percentiles[i] && is_dual_plane && permit_encode)
				bsd->decimation_mode_percentile_2planes[decimation_mode] = percentiles[i];
		}

	}

	if (xdim * ydim * zdim <= 64)
	{
		bsd->texelcount_for_bitmap_partitioning = xdim * ydim * zdim;
		for (i = 0; i < xdim * ydim * zdim; i++)
			bsd->texels_for_bitmap_partitioning[i] = i;
	}

	else
	{
		// pick 64 random texels for use with bitmap partitioning.
		int arr[MAX_TEXELS_PER_BLOCK];
		for (i = 0; i < xdim * ydim * zdim; i++)
			arr[i] = 0;
		int arr_elements_set = 0;
		while (arr_elements_set < 64)
		{
			int idx = rand() % (xdim * ydim * zdim);
			if (arr[idx] == 0)
			{
				arr_elements_set++;
				arr[idx] = 1;
			}
		}
		int texel_weights_written = 0;
		int idx = 0;
		while (texel_weights_written < 64)
		{
			if (arr[idx])
				bsd->texels_for_bitmap_partitioning[texel_weights_written++] = idx;
			idx++;
		}
		bsd->texelcount_for_bitmap_partitioning = 64;
	}
}


static void sort_decimation_modes(const block_size_descriptor * bsd, int is_dual_plane, block_size_descriptor_sorted *sorted_bsd, int * decimation_mode_mapping)
{
	float percentiles[MAX_DECIMATION_MODES];
	int valid_modes = 0;

	for (int i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		percentiles[i] = is_dual_plane ? bsd->decimation_mode_percentile_2planes[i] : bsd->decimation_mode_percentile_1plane[i];

		if (percentiles[i] <= 1.0f)
			valid_modes++;
	}

	for (int i = 0; i < valid_modes; i++)
	{
		float lowest_percentile = INFINITY;
		int best_mode = -1;
		for (int j = 0; j < MAX_DECIMATION_MODES; j++)
		{
			if (percentiles[j] < lowest_percentile)
			{
				lowest_percentile = percentiles[j];
				best_mode = j;
			}
		}

		percentiles[best_mode] = INFINITY;
		decimation_mode_mapping[best_mode] = i;
		int maxprec = is_dual_plane ? bsd->decimation_mode_maxprec_2planes[best_mode] : bsd->decimation_mode_maxprec_1plane[best_mode];
		sorted_bsd->decimation_mode_maxprec[i] = maxprec;
		sorted_bsd->decimation_mode_percentile[i] = lowest_percentile;
		sorted_bsd->decimation_mode_samples[i] = bsd->decimation_mode_samples[best_mode];
		sorted_bsd->decimation_tables[i] = bsd->decimation_tables[best_mode];
	}
}

// find weight mode with lowest percentile
static int find_best_weight_mode(const float *percentiles, const block_size_descriptor *bsd, int is_dual_plane)
{
	float lowest_percentile = INFINITY;
	int best_weight_mode = -1;

	for (int j = 0; j < MAX_WEIGHT_MODES; j++)
	{
		if ((bsd->block_modes[j].is_dual_plane != is_dual_plane) || !bsd->block_modes[j].permit_encode)
			continue;

		if (percentiles[j] < lowest_percentile)
		{
			lowest_percentile = percentiles[j];
			best_weight_mode = j;
		}
	}

	return best_weight_mode;
}

static void sort_weight_modes(const block_size_descriptor * bsd, int is_dual_plane, block_size_descriptor_sorted *sorted_bsd, const int * decimation_mode_mapping)
{
	float percentiles[MAX_WEIGHT_MODES];
	int valid_modes = 0;

	for (int i = 0; i < MAX_WEIGHT_MODES; i++)
	{
		percentiles[i] = bsd->block_modes[i].percentile;

		if ((bsd->block_modes[i].permit_encode) && (is_dual_plane == bsd->block_modes[i].is_dual_plane))
			valid_modes++;
	}

	assert(valid_modes >= 4);

	int best_weight_mode = find_best_weight_mode(percentiles, bsd, is_dual_plane);
	int best_decimation_mode = bsd->block_modes[best_weight_mode].decimation_mode;
	int num_weights = bsd->decimation_tables[best_decimation_mode]->num_weights;
	if (is_dual_plane)
		num_weights *= 2;
	// max available bits:
	// 76 in 3partition 2planes mode
	// 70 in 4partition 1plane mode
	const int max_weight_bits = (is_dual_plane) ? 76 : 70;

	for (int i = 0; i < 4; i++)
	{
		int best_block_mode = -1;
		float lowest_percentile = INFINITY;

		for (int j = 0; j < MAX_WEIGHT_MODES; j++)
		{
			if ((bsd->block_modes[j].is_dual_plane != is_dual_plane)
				|| (bsd->block_modes[j].decimation_mode != best_decimation_mode)
				|| !bsd->block_modes[j].permit_encode)
				continue;

			int bits_used_by_weights = compute_ise_bitcount(num_weights, (quantization_method)bsd->block_modes[j].quantization_mode);

			if (bits_used_by_weights > max_weight_bits)
				continue;

			if (percentiles[j] < lowest_percentile)
			{
				lowest_percentile = percentiles[j];
				best_block_mode = j;
			}
		}

		percentiles[best_block_mode] = INFINITY;
		sorted_bsd->block_modes[i].block_mode = best_block_mode;
		sorted_bsd->block_modes[i].quantization_mode = bsd->block_modes[best_block_mode].quantization_mode;
		sorted_bsd->block_modes[i].sorted_decimation_mode = decimation_mode_mapping[best_decimation_mode];
	}

	for (int i = 4; i < valid_modes; i++)
	{
		int best_block_mode = find_best_weight_mode(percentiles, bsd, is_dual_plane);

		percentiles[best_block_mode] = INFINITY;
		int old_decimation_mode = bsd->block_modes[best_block_mode].decimation_mode;
		sorted_bsd->block_modes[i].block_mode = best_block_mode;
		sorted_bsd->block_modes[i].quantization_mode = bsd->block_modes[best_block_mode].quantization_mode;
		sorted_bsd->block_modes[i].sorted_decimation_mode = decimation_mode_mapping[old_decimation_mode];
	}
}

int get_decimation_mode_limit(int xdim, int ydim, int zdim, int is_dual_plane, float mode_cutoff)
{
	const block_size_descriptor * bsd = get_block_size_descriptor(xdim, ydim, zdim);
	const block_size_descriptor_sorted * sorted_bsd = get_sorted_block_size_descriptor(xdim, ydim, zdim, is_dual_plane);

	int i;

	for (i = 0; i < MAX_DECIMATION_MODES; i++)
		if (sorted_bsd->decimation_mode_maxprec[i] < 0 || sorted_bsd->decimation_mode_percentile[i] > mode_cutoff)
			break;

	assert(i > 0);

	return i;
}

int get_weight_mode_limit(int xdim, int ydim, int zdim, int is_dual_plane, float mode_cutoff)
{
	const block_size_descriptor * bsd = get_block_size_descriptor(xdim, ydim, zdim);
	const block_size_descriptor_sorted * sorted_bsd = get_sorted_block_size_descriptor(xdim, ydim, zdim, is_dual_plane);

	int i;

	for (i = 0; i < MAX_SORTED_WEIGHT_MODES; i++)
	{
		int block_mode = sorted_bsd->block_modes[i].block_mode;
		if (block_mode < 0)
			break;
		if (i > 3 && bsd->block_modes[block_mode].percentile > mode_cutoff)
			break;
	}

	// during compression we always test at least 4 modes
	assert(i >= 4);

	return i;
}

static void construct_sorted_block_size_descriptors(const block_size_descriptor * bsd, block_size_descriptor_sorted *sorted_bsd, int is_dual_plane)
{
	int decimation_mode_map[MAX_DECIMATION_MODES];

	memset(sorted_bsd, -1, sizeof(block_size_descriptor_sorted));

	sort_decimation_modes(bsd, is_dual_plane, sorted_bsd, decimation_mode_map);
	sort_weight_modes(bsd, is_dual_plane, sorted_bsd, decimation_mode_map);
}


// function to obtain a block size descriptor. If the descriptor does not exist,
// it is created as needed. Should not be called from within multithreaded code.
const block_size_descriptor *get_block_size_descriptor(int xdim, int ydim, int zdim)
{
	int bsd_index = xdim + (ydim << 4) + (zdim << 8);
	if (bsd_pointers[bsd_index] == NULL)
	{
		block_size_info *bsd_info = new block_size_info;
		block_size_descriptor *bsd = &bsd_info->bsd;
		if (zdim > 1)
			construct_block_size_descriptor_3d(xdim, ydim, zdim, bsd);
		else
			construct_block_size_descriptor_2d(xdim, ydim, bsd);

		construct_sorted_block_size_descriptors(bsd, &bsd_info->bsd_1plane, 0);
		construct_sorted_block_size_descriptors(bsd, &bsd_info->bsd_2planes, 1);

		bsd_pointers[bsd_index] = bsd_info;
	}
	return &bsd_pointers[bsd_index]->bsd;
}

const block_size_descriptor_sorted *get_sorted_block_size_descriptor(int xdim, int ydim, int zdim, int is_dual_plane)
{
	int bsd_index = xdim + (ydim << 4) + (zdim << 8);
	if (bsd_pointers[bsd_index] == NULL)
	{
		// construct bsd
		get_block_size_descriptor(xdim, ydim, zdim);
	}

	block_size_info *bsd_info = bsd_pointers[bsd_index];

	return is_dual_plane ? &(bsd_info->bsd_2planes) : &(bsd_info->bsd_1plane);
}
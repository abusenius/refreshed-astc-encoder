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
 *	@brief	Compress a block of colors, expressed as a symbolic block, for ASTC.
 */ 
/*----------------------------------------------------------------------------*/ 

#include "astc_codec_internals.h"
#include "astc_codec_batch.h"

#include "softfloat.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifdef DEBUG_CAPTURE_NAN
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE
	#endif

	#include <fenv.h>
#endif

#include <stdio.h>

int realign_weights(astc_decode_mode decode_mode,
					int xdim, int ydim, int zdim, const imageblock * blk, const error_weight_block * ewb, symbolic_compressed_block * scb, uint8_t * weight_set8, uint8_t * plane2_weight_set8)
{
	int i, j;

	// get the appropriate partition descriptor.
	int partition_count = scb->partition_count;
	const partition_info *pt = get_partition_table(xdim, ydim, zdim, partition_count);
	pt += scb->partition_index;

	// get the appropriate block descriptor
	const block_size_descriptor *bsd = get_block_size_descriptor(xdim, ydim, zdim);
	const decimation_table *const *ixtab2 = bsd->decimation_tables;

	const decimation_table *it = ixtab2[bsd->block_modes[scb->block_mode].decimation_mode];

	int is_dual_plane = bsd->block_modes[scb->block_mode].is_dual_plane;

	// get quantization-parameters
	int weight_quantization_level = bsd->block_modes[scb->block_mode].quantization_mode;


	// decode the color endpoints
	ushort4 color_endpoint0[4];
	ushort4 color_endpoint1[4];
	int rgb_hdr[4];
	int alpha_hdr[4];
	int nan_endpoint[4];


	for (i = 0; i < partition_count; i++)
		unpack_color_endpoints(decode_mode,
							   scb->color_formats[i], scb->color_quantization_level, scb->color_values[i], &rgb_hdr[i], &alpha_hdr[i], &nan_endpoint[i], &(color_endpoint0[i]), &(color_endpoint1[i]));


	float uq_plane1_weights[MAX_WEIGHTS_PER_BLOCK];
	float uq_plane2_weights[MAX_WEIGHTS_PER_BLOCK];
	int weight_count = it->num_weights;

	// read the weights.

	const quantization_and_transfer_table *qat = &(quant_and_xfer_tables[weight_quantization_level]);

	for (i = 0; i < weight_count; i++)
	{
		uq_plane1_weights[i] = ((float)weight_set8[i]) / 64;
	}
	if (is_dual_plane)
	{
		for (i = 0; i < weight_count; i++)
			uq_plane2_weights[i] = ((float)plane2_weight_set8[i]) / 64;
	}


	int plane2_color_component = is_dual_plane ? scb->plane2_color_component : -1;

	// for each weight, unquantize the weight, use it to compute a color and a color error.
	// then, increment the weight until the color error stops decreasing
	// then, decrement the weight until the color error stops increasing
	
	#define COMPUTE_ERROR( errorvar ) \
		errorvar = 0.0f; \
		for(j=0;j<texels_to_evaluate;j++) \
			{ \
			int texel = it->weight_texel[i][j]; \
			int partition = pt->partition_of_texel[texel]; \
			float plane1_weight = compute_value_of_texel_flt( texel, it, uq_plane1_weights ); \
			float plane2_weight = 0.0f; \
			if( is_dual_plane ) \
				plane2_weight = compute_value_of_texel_flt( texel, it, uq_plane2_weights ); \
			int int_plane1_weight = static_cast<int>(floor( plane1_weight*64.0f + 0.5f ) ); \
			int int_plane2_weight = static_cast<int>(floor( plane2_weight*64.0f + 0.5f ) ); \
			ushort4 lrp_color = lerp_color_int( \
				decode_mode, \
				color_endpoint0[partition], \
				color_endpoint1[partition], \
				int_plane1_weight, \
				int_plane2_weight, \
				plane2_color_component ); \
			float4 color = float4( lrp_color.x, lrp_color.y, lrp_color.z, lrp_color.w ); \
			float4 origcolor = float4( \
				blk->work_data[4*texel], \
				blk->work_data[4*texel+1], \
				blk->work_data[4*texel+2], \
				blk->work_data[4*texel+3] ); \
			float4 error_weight = ewb->error_weights[texel]; \
			float4 colordiff = color - origcolor; \
			errorvar += dot( colordiff*colordiff, error_weight ); \
			}


	int adjustments = 0;

	for (i = 0; i < weight_count; i++)
	{
		int current_wt = weight_set8[i];
		int texels_to_evaluate = it->weight_num_texels[i];

		float current_error;

		COMPUTE_ERROR(current_error);

		// increment until error starts increasing.
		while (1)
		{
			if (64 == current_wt)
				break;
			int next_wt = qat->next_unquantized_value[current_wt];
			uq_plane1_weights[i] = ((float)next_wt) / 64;
			float next_error;
			COMPUTE_ERROR(next_error);
			if (next_error < current_error)
			{
				// succeeded, increment the weight
				current_wt = next_wt;
				current_error = next_error;
				adjustments++;
			}
			else
			{
				// failed, back out the attempted increment
				uq_plane1_weights[i] = ((float)current_wt) / 64;
				break;
			}
		}
		// decrement until error starts increasing
		while (1)
		{
			if (0 == current_wt)
				break;
			int prev_wt = qat->prev_unquantized_value[current_wt];
			uq_plane1_weights[i] = ((float)prev_wt) / 64;
			float prev_error;
			COMPUTE_ERROR(prev_error);
			if (prev_error < current_error)
			{
				// succeeded, decrement the weight
				current_wt = prev_wt;
				current_error = prev_error;
				adjustments++;
			}
			else
			{
				// failed, back out the attempted decrement
				uq_plane1_weights[i] = ((float)current_wt) / 64;
				break;
			}
		}

		weight_set8[i] = current_wt;
	}

	if (!is_dual_plane)
		return adjustments;

	// processing of the second plane of weights
	for (i = 0; i < weight_count; i++)
	{
		int current_wt = plane2_weight_set8[i];
		int texels_to_evaluate = it->weight_num_texels[i];

		float current_error;

		COMPUTE_ERROR(current_error);

		// increment until error starts increasing.
		while (1)
		{
			if (64 == current_wt)
				break;
			int next_wt = qat->next_unquantized_value[current_wt];
			uq_plane2_weights[i] = ((float)next_wt) / 64;
			float next_error;
			COMPUTE_ERROR(next_error);
			if (next_error < current_error)
			{
				// succeeded, increment the weight
				current_wt = next_wt;
				current_error = next_error;
				adjustments++;
			}
			else
			{
				// failed, back out the attempted increment
				uq_plane2_weights[i] = ((float)current_wt) / 64;
				break;
			}
		}
		// decrement until error starts increasing
		while (1)
		{
			if (0 == current_wt)
				break;
			int prev_wt = qat->prev_unquantized_value[current_wt];
			uq_plane1_weights[i] = ((float)prev_wt) / 64;
			float prev_error;
			COMPUTE_ERROR(prev_error);
			if (prev_error < current_error)
			{
				// succeeded, decrement the weight
				current_wt = prev_wt;
				current_error = prev_error;
				adjustments++;
			}
			else
			{
				// failed, back out the attempted decrement
				uq_plane2_weights[i] = ((float)current_wt) / 64;
				break;
			}
		}

		plane2_weight_set8[i] = current_wt;
	}

	return adjustments;
}




void expand_block_artifact_suppression(int xdim, int ydim, int zdim, error_weighting_params * ewp)
{
	int x, y, z;
	float centerpos_x = (xdim - 1) * 0.5f;
	float centerpos_y = (ydim - 1) * 0.5f;
	float centerpos_z = (zdim - 1) * 0.5f;
	float *bef = ewp->block_artifact_suppression_expanded;

	for (z = 0; z < zdim; z++)
		for (y = 0; y < ydim; y++)
			for (x = 0; x < xdim; x++)
			{
				float xdif = (x - centerpos_x) / xdim;
				float ydif = (y - centerpos_y) / ydim;
				float zdif = (z - centerpos_z) / zdim;

				float wdif = 0.36f;
				float dist = sqrt(xdif * xdif + ydif * ydif + zdif * zdif + wdif * wdif);
				*bef = pow(dist, ewp->block_artifact_suppression);
				bef++;
			}
}



// Function to set error weights for each color component for each texel in a block.
// Returns the sum of all the error values set.

float prepare_error_weight_block(const astc_codec_image * input_image,
								 int xdim, int ydim, int zdim, const error_weighting_params * ewp, const imageblock * blk, error_weight_block * ewb)
{

	int x, y, z;
	int idx = 0;

	int any_mean_stdev_weight =
		ewp->rgb_base_weight != 1.0 || ewp->alpha_base_weight != 1.0 || ewp->rgb_mean_weight != 0.0 || ewp->rgb_stdev_weight != 0.0 || ewp->alpha_mean_weight != 0.0 || ewp->alpha_stdev_weight != 0.0;

	float4 color_weights = float4(ewp->rgba_weights[0],
								  ewp->rgba_weights[1],
								  ewp->rgba_weights[2],
								  ewp->rgba_weights[3]);

	ewb->contains_zeroweight_texels = 0;

	for (z = 0; z < zdim; z++)
		for (y = 0; y < ydim; y++)
			for (x = 0; x < xdim; x++)
			{
				int xpos = x + blk->xpos;
				int ypos = y + blk->ypos;
				int zpos = z + blk->zpos;

				if (xpos >= input_image->xsize || ypos >= input_image->ysize || zpos >= input_image->zsize)
				{
					float4 weights = float4(1e-11f, 1e-11f, 1e-11f, 1e-11f);
					ewb->error_weights[idx] = weights;
					ewb->contains_zeroweight_texels = 1;
				}
				else
				{
					float4 error_weight = float4(ewp->rgb_base_weight,
												 ewp->rgb_base_weight,
												 ewp->rgb_base_weight,
												 ewp->alpha_base_weight);

					if (any_mean_stdev_weight)
					{
						float4 avg = input_averages[zpos][ypos][xpos];
						if (avg.x < 6e-5f)
							avg.x = 6e-5f;
						if (avg.y < 6e-5f)
							avg.y = 6e-5f;
						if (avg.z < 6e-5f)
							avg.z = 6e-5f;
						if (avg.w < 6e-5f)
							avg.w = 6e-5f;
						/* 
						   printf("avg: %f %f %f %f\n", avg.x, avg.y, avg.z, avg.w ); */
						avg = avg * avg;

						float4 variance = input_variances[zpos][ypos][xpos];
						variance = variance * variance;

						float favg = (avg.x + avg.y + avg.z) * (1.0f / 3.0f);
						float fvar = (variance.x + variance.y + variance.z) * (1.0f / 3.0f);

						float mixing = ewp->rgb_mean_and_stdev_mixing;
						avg.xyz = float3(favg, favg, favg) * mixing + avg.xyz * (1.0f - mixing);
						variance.xyz = float3(fvar, fvar, fvar) * mixing + variance.xyz * (1.0f - mixing);

						float4 stdev = float4(sqrt(MAX(variance.x, 0.0f)),
											  sqrt(MAX(variance.y, 0.0f)),
											  sqrt(MAX(variance.z, 0.0f)),
											  sqrt(MAX(variance.w, 0.0f)));

						avg.xyz = avg.xyz * ewp->rgb_mean_weight;
						avg.w = avg.w * ewp->alpha_mean_weight;
						stdev.xyz = stdev.xyz * ewp->rgb_stdev_weight;
						stdev.w = stdev.w * ewp->alpha_stdev_weight;
						error_weight = error_weight + avg + stdev;

						error_weight = float4(1.0f, 1.0f, 1.0f, 1.0f) / error_weight;
					}

					if (ewp->ra_normal_angular_scale)
					{
						float x = (blk->orig_data[4 * idx] - 0.5f) * 2.0f;
						float y = (blk->orig_data[4 * idx + 3] - 0.5f) * 2.0f;
						float denom = 1.0f - x * x - y * y;
						if (denom < 0.1f)
							denom = 0.1f;
						denom = 1.0f / denom;
						error_weight.x *= 1.0f + x * x * denom;
						error_weight.w *= 1.0f + y * y * denom;
					}

					if (ewp->enable_rgb_scale_with_alpha)
					{
						float alpha_scale;
						if (ewp->alpha_radius != 0)
							alpha_scale = input_alpha_averages[zpos][ypos][xpos];
						else
							alpha_scale = blk->orig_data[4 * idx + 3];
						if (alpha_scale < 0.0001f)
							alpha_scale = 0.0001f;
						alpha_scale *= alpha_scale;
						error_weight.xyz = error_weight.xyz * alpha_scale;
					}
					error_weight = error_weight * color_weights;
					error_weight = error_weight * ewp->block_artifact_suppression_expanded[idx];

					// if we perform a conversion from linear to sRGB, then we multiply
					// the weight with the derivative of the linear->sRGB transform function.
					if (perform_srgb_transform)
					{
						float r = blk->orig_data[4 * idx];
						float g = blk->orig_data[4 * idx + 1];
						float b = blk->orig_data[4 * idx + 2];
						if (r < 0.0031308f)
							r = 12.92f;
						else
							r = 0.4396f * pow(r, -0.58333f);
						if (g < 0.0031308f)
							g = 12.92f;
						else
							g = 0.4396f * pow(g, -0.58333f);
						if (b < 0.0031308f)
							b = 12.92f;
						else
							b = 0.4396f * pow(b, -0.58333f);
						error_weight.x *= r;
						error_weight.y *= g;
						error_weight.z *= b;
					}

					/*
						printf("%f %f %f %f\n", error_weight.x, error_weight.y, error_weight.z, error_weight.w );
					*/

					// when we loaded the block to begin with, we applied a transfer function
					// and computed the derivative of the transfer function. However, the
					// error-weight computation so far is based on the original color values,
					// not the transfer-function values. As such, we must multiply the
					// error weights by the derivative of the inverse of the transfer function,
					// which is equivalent to dividing by the derivative of the transfer
					// function.

					error_weight.x /= (blk->deriv_data[4 * idx] * blk->deriv_data[4 * idx] * 1e-10f);
					error_weight.y /= (blk->deriv_data[4 * idx + 1] * blk->deriv_data[4 * idx + 1] * 1e-10f);
					error_weight.z /= (blk->deriv_data[4 * idx + 2] * blk->deriv_data[4 * idx + 2] * 1e-10f);
					error_weight.w /= (blk->deriv_data[4 * idx + 3] * blk->deriv_data[4 * idx + 3] * 1e-10f);

					/*
						printf("--> %f %f %f %f\n", error_weight.x, error_weight.y, error_weight.z, error_weight.w );
					*/

					ewb->error_weights[idx] = error_weight;
					if (dot(error_weight, float4(1, 1, 1, 1)) < 1e-10f)
						ewb->contains_zeroweight_texels = 1;
				}
				idx++;
			}

	int i;

	float4 error_weight_sum = float4(0, 0, 0, 0);
	int texels_per_block = xdim * ydim * zdim;

	for (i = 0; i < texels_per_block; i++)
	{
		error_weight_sum = error_weight_sum + ewb->error_weights[i];

		ewb->texel_weight_r[i] = ewb->error_weights[i].x;
		ewb->texel_weight_g[i] = ewb->error_weights[i].y;
		ewb->texel_weight_b[i] = ewb->error_weights[i].z;
		ewb->texel_weight_a[i] = ewb->error_weights[i].w;

		ewb->texel_weight_rg[i] = (ewb->error_weights[i].x + ewb->error_weights[i].y) * 0.5f;
		ewb->texel_weight_rb[i] = (ewb->error_weights[i].x + ewb->error_weights[i].z) * 0.5f;
		ewb->texel_weight_gb[i] = (ewb->error_weights[i].y + ewb->error_weights[i].z) * 0.5f;
		ewb->texel_weight_ra[i] = (ewb->error_weights[i].x + ewb->error_weights[i].w) * 0.5f;

		ewb->texel_weight_gba[i] = (ewb->error_weights[i].y + ewb->error_weights[i].z + ewb->error_weights[i].w) * 0.333333f;
		ewb->texel_weight_rba[i] = (ewb->error_weights[i].x + ewb->error_weights[i].z + ewb->error_weights[i].w) * 0.333333f;
		ewb->texel_weight_rga[i] = (ewb->error_weights[i].x + ewb->error_weights[i].y + ewb->error_weights[i].w) * 0.333333f;
		ewb->texel_weight_rgb[i] = (ewb->error_weights[i].x + ewb->error_weights[i].y + ewb->error_weights[i].z) * 0.333333f;
		ewb->texel_weight[i] = (ewb->error_weights[i].x + ewb->error_weights[i].y + ewb->error_weights[i].z + ewb->error_weights[i].w) * 0.25f;
	}

	return dot(error_weight_sum, float4(1, 1, 1, 1));
}


/* 
	functions to analyze block statistical properties:
		* simple properties: * mean * variance
		* covariance-matrix correllation coefficients
 */


// compute averages and covariance matrices for 4 components
static void compute_covariance_matrix(int xdim, int ydim, int zdim, const imageblock * blk, const error_weight_block * ewb, mat4 * cov_matrix)
{
	int i;

	int texels_per_block = xdim * ydim * zdim;

	float r_sum = 0.0f;
	float g_sum = 0.0f;
	float b_sum = 0.0f;
	float a_sum = 0.0f;
	float rr_sum = 0.0f;
	float gg_sum = 0.0f;
	float bb_sum = 0.0f;
	float aa_sum = 0.0f;
	float rg_sum = 0.0f;
	float rb_sum = 0.0f;
	float ra_sum = 0.0f;
	float gb_sum = 0.0f;
	float ga_sum = 0.0f;
	float ba_sum = 0.0f;

	float weight_sum = 0.0f;

	for (i = 0; i < texels_per_block; i++)
	{
		float weight = ewb->texel_weight[i];
		if (weight < 0.0f)
			ASTC_CODEC_INTERNAL_ERROR;
		weight_sum += weight;
		float r = blk->work_data[4 * i];
		float g = blk->work_data[4 * i + 1];
		float b = blk->work_data[4 * i + 2];
		float a = blk->work_data[4 * i + 3];
		r_sum += r * weight;
		rr_sum += r * (r * weight);
		rg_sum += g * (r * weight);
		rb_sum += b * (r * weight);
		ra_sum += a * (r * weight);
		g_sum += g * weight;
		gg_sum += g * (g * weight);
		gb_sum += b * (g * weight);
		ga_sum += a * (g * weight);
		b_sum += b * weight;
		bb_sum += b * (b * weight);
		ba_sum += a * (b * weight);
		a_sum += a * weight;
		aa_sum += a * (a * weight);
	}

	float rpt = 1.0f / MAX(weight_sum, 1e-7f);
	float rs = r_sum;
	float gs = g_sum;
	float bs = b_sum;
	float as = a_sum;

	cov_matrix->v[0] = float4(rr_sum - rs * rs * rpt, rg_sum - rs * gs * rpt, rb_sum - rs * bs * rpt, ra_sum - rs * as * rpt);
	cov_matrix->v[1] = float4(rg_sum - rs * gs * rpt, gg_sum - gs * gs * rpt, gb_sum - gs * bs * rpt, ga_sum - gs * as * rpt);
	cov_matrix->v[2] = float4(rb_sum - rs * bs * rpt, gb_sum - gs * bs * rpt, bb_sum - bs * bs * rpt, ba_sum - bs * as * rpt);
	cov_matrix->v[3] = float4(ra_sum - rs * as * rpt, ga_sum - gs * as * rpt, ba_sum - bs * as * rpt, aa_sum - as * as * rpt);

}



void prepare_block_statistics(int xdim, int ydim, int zdim, const imageblock * blk, const error_weight_block * ewb, int *is_normal_map, float *lowest_correl)
{
	int i;

	mat4 cov_matrix;

	compute_covariance_matrix(xdim, ydim, zdim, blk, ewb, &cov_matrix);

	// use the covariance matrix to compute
	// correllation coefficients
	float rr_var = cov_matrix.v[0].x;
	float gg_var = cov_matrix.v[1].y;
	float bb_var = cov_matrix.v[2].z;
	float aa_var = cov_matrix.v[3].w;

	float rg_correlation = cov_matrix.v[0].y / sqrt(MAX(rr_var * gg_var, 1e-30f));
	float rb_correlation = cov_matrix.v[0].z / sqrt(MAX(rr_var * bb_var, 1e-30f));
	float ra_correlation = cov_matrix.v[0].w / sqrt(MAX(rr_var * aa_var, 1e-30f));
	float gb_correlation = cov_matrix.v[1].z / sqrt(MAX(gg_var * bb_var, 1e-30f));
	float ga_correlation = cov_matrix.v[1].w / sqrt(MAX(gg_var * aa_var, 1e-30f));
	float ba_correlation = cov_matrix.v[2].w / sqrt(MAX(bb_var * aa_var, 1e-30f));

	if (astc_isnan(rg_correlation))
		rg_correlation = 1.0f;
	if (astc_isnan(rb_correlation))
		rb_correlation = 1.0f;
	if (astc_isnan(ra_correlation))
		ra_correlation = 1.0f;
	if (astc_isnan(gb_correlation))
		gb_correlation = 1.0f;
	if (astc_isnan(ga_correlation))
		ga_correlation = 1.0f;
	if (astc_isnan(ba_correlation))
		ba_correlation = 1.0f;

	float lowest_correlation = MIN(fabs(rg_correlation), fabs(rb_correlation));
	lowest_correlation = MIN(lowest_correlation, fabs(ra_correlation));
	lowest_correlation = MIN(lowest_correlation, fabs(gb_correlation));
	lowest_correlation = MIN(lowest_correlation, fabs(ga_correlation));
	lowest_correlation = MIN(lowest_correlation, fabs(ba_correlation));
	*lowest_correl = lowest_correlation;

	// compute a "normal-map" factor
	// this factor should be exactly 0.0 for a normal map, while it may be all over the
	// place for anything that is NOT a normal map. We can probably assume that a factor
	// of less than 0.2f represents a normal map.

	float nf_sum = 0.0f;

	int texels_per_block = xdim * ydim * zdim;

	for (i = 0; i < texels_per_block; i++)
	{
		float3 val = float3(blk->orig_data[4 * i],
							blk->orig_data[4 * i + 1],
							blk->orig_data[4 * i + 2]);
		val = (val - float3(0.5f, 0.5f, 0.5f)) * 2.0f;
		float length_squared = dot(val, val);
		float nf = fabs(length_squared - 1.0f);
		nf_sum += nf;
	}
	float nf_avg = nf_sum / texels_per_block;
	*is_normal_map = nf_avg < 0.2;
}





void compress_constant_color_block(int xdim, int ydim, int zdim, const imageblock * blk, const error_weight_block * ewb, symbolic_compressed_block * scb)
{
	int texel_count = xdim * ydim * zdim;
	int i;

	float4 color_sum = float4(0, 0, 0, 0);
	float4 color_weight_sum = float4(0, 0, 0, 0);

	const float *clp = blk->work_data;
	for (i = 0; i < texel_count; i++)
	{
		float4 weights = ewb->error_weights[i];
		float4 color_data = float4(clp[4 * i], clp[4 * i + 1], clp[4 * i + 2], clp[4 * i + 3]);
		color_sum = color_sum + (color_data * weights);
		color_weight_sum = color_weight_sum + weights;
	}

	float4 avg_color = color_sum / color_weight_sum;

	int use_fp16 = blk->rgb_lns[0];

	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
		{
			printf("Averaged color: %f %f %f %f\n", avg_color.x, avg_color.y, avg_color.z, avg_color.w);
		}
	#endif

	// convert the color
	if (blk->rgb_lns[0])
	{
		int avg_red = static_cast < int >(floor(avg_color.x + 0.5f));
		int avg_green = static_cast < int >(floor(avg_color.y + 0.5f));
		int avg_blue = static_cast < int >(floor(avg_color.z + 0.5f));

		if (avg_red < 0)
			avg_red = 0;
		else if (avg_red > 65535)
			avg_red = 65535;

		if (avg_green < 0)
			avg_green = 0;
		else if (avg_green > 65535)
			avg_green = 65535;

		if (avg_blue < 0)
			avg_blue = 0;
		else if (avg_blue > 65535)
			avg_blue = 65535;

		avg_color.x = sf16_to_float(lns_to_sf16(avg_red));
		avg_color.y = sf16_to_float(lns_to_sf16(avg_green));
		avg_color.z = sf16_to_float(lns_to_sf16(avg_blue));
	}
	else
	{
		avg_color.x *= (1.0f / 65535.0f);
		avg_color.y *= (1.0f / 65535.0f);
		avg_color.z *= (1.0f / 65535.0f);
	}
	if (blk->alpha_lns[0])
	{
		int avg_alpha = static_cast < int >(floor(avg_color.w + 0.5f));

		if (avg_alpha < 0)
			avg_alpha = 0;
		else if (avg_alpha > 65535)
			avg_alpha = 65535;

		avg_color.w = sf16_to_float(lns_to_sf16(avg_alpha));
	}
	else
	{
		avg_color.w *= (1.0f / 65535.0f);
	}

#ifdef DEBUG_PRINT_DIAGNOSTICS
	if (print_diagnostics)
	{
		printf("Averaged color: %f %f %f %f   (%d)\n", avg_color.x, avg_color.y, avg_color.z, avg_color.w, use_fp16);

	}
#endif

	if (use_fp16)
	{
		scb->error_block = 0;
		scb->block_mode = -1;
		scb->partition_count = 0;
		scb->constant_color[0] = float_to_sf16(avg_color.x, SF_NEARESTEVEN);
		scb->constant_color[1] = float_to_sf16(avg_color.y, SF_NEARESTEVEN);
		scb->constant_color[2] = float_to_sf16(avg_color.z, SF_NEARESTEVEN);
		scb->constant_color[3] = float_to_sf16(avg_color.w, SF_NEARESTEVEN);
	}

	else
	{
		scb->error_block = 0;
		scb->block_mode = -2;
		scb->partition_count = 0;
		float red = avg_color.x;
		float green = avg_color.y;
		float blue = avg_color.z;
		float alpha = avg_color.w;
		if (red < 0)
			red = 0;
		else if (red > 1)
			red = 1;
		if (green < 0)
			green = 0;
		else if (green > 1)
			green = 1;
		if (blue < 0)
			blue = 0;
		else if (blue > 1)
			blue = 1;
		if (alpha < 0)
			alpha = 0;
		else if (alpha > 1)
			alpha = 1;
		scb->constant_color[0] = static_cast < int >(floor(red * 65535.0f + 0.5f));
		scb->constant_color[1] = static_cast < int >(floor(green * 65535.0f + 0.5f));
		scb->constant_color[2] = static_cast < int >(floor(blue * 65535.0f + 0.5f));
		scb->constant_color[3] = static_cast < int >(floor(alpha * 65535.0f + 0.5f));
	}
}

//SymbolicBatchCompressor::SymbolicBatchCompressor(int _max_batch_size)
//	: max_batch_size(_max_batch_size), xdim(4), ydim(4), zdim(1), decode_mode(astc_decode_mode::DECODE_LDR)
//{
//	allocate_buffers(max_batch_size);
//}

SymbolicBatchCompressor::SymbolicBatchCompressor(int _max_batch_size, int _xdim, int _ydim, int _zdim, astc_decode_mode _decode_mode, const error_weighting_params * _ewp)
	: max_batch_size(_max_batch_size), xdim(_xdim), ydim(_ydim), zdim(_zdim), decode_mode(_decode_mode)
{
	ewp = *_ewp;

	cl_int status;
	opencl_queue = clCreateCommandQueue(opencl_context, opencl_device, 0, &status);
	OCL_CHECK_STATUS("Unable to create command queue");

	const partition_statistics * pstat;
	
	for (size_t pcount = 2; pcount <= 4; pcount++)
	{
		pstat = get_partition_stats(xdim, ydim, zdim, pcount);
		fbp.partition_search_limits[pcount] = MIN(ewp.partition_search_limit, pstat->unique_partitionings_with_all_partitions);
	}

	int texels_per_block = xdim * ydim * zdim;
	fbp.weight_imprecision_estim_squared = calculate_weight_imprecision(texels_per_block);

	allocate_buffers(max_batch_size);

	OCL_CREATE_BUFFER(blk_buf, CL_MEM_READ_ONLY, sizeof(imageblock) * max_batch_size, NULL);
	
	for (size_t pcount = 2; pcount <= 4; pcount++)
	{
		OCL_CREATE_BUFFER(fbp.ptab[pcount], CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, sizeof(partition_info) * PARTITION_COUNT, NULL);
		OCL_WRITE_BUFFER(fbp.ptab[pcount], sizeof(partition_info) * PARTITION_COUNT, get_partition_table(xdim, ydim, zdim, pcount));
	}

	OCL_CREATE_KERNEL(fbp, find_best_partitionings);
	OCL_SET_KERNEL_ARG(fbp.find_best_partitionings, 0, blk_stat);
	OCL_SET_KERNEL_ARG(fbp.find_best_partitionings, 1, blk_buf);
	OCL_SET_KERNEL_ARG(fbp.find_best_partitionings, 2, fbp.partition_sequence);
	OCL_SET_KERNEL_ARG(fbp.find_best_partitionings, 3, partition_indices_1plane_batch);
	OCL_SET_KERNEL_ARG(fbp.find_best_partitionings, 4, partition_indices_2planes_batch);
	OCL_SET_KERNEL_ARG(fbp.find_best_partitionings, 5, ewb_batch);
	
	clFlush(opencl_queue);
}

#define FIND_BEST_SCB_CANDIDATES(modesel) {best_errorval_in_mode = 1e30f;\
		for (j = 0; j < SCB_CANDIDATES; j++)\
		{\
			imageblock temp;\
			auto scb2test = &scb_candidates[SCB_CANDIDATES * blk_idx];\
			if (scb2test[j].error_block)\
				continue;\
			decompress_symbolic_block(decode_mode, xdim, ydim, zdim, 0, 0, 0, scb2test + j, &temp);\
			float errorval = compute_imageblock_difference(xdim, ydim, zdim, &blk_batch[blk_idx], &temp, &ewb_batch[blk_idx]);\
		\
			if (errorval < best_errorval_in_mode)\
				best_errorval_in_mode = errorval;\
		\
			if (errorval < error_of_best_block[blk_idx])\
			{\
				error_of_best_block[blk_idx] = errorval;\
				scb_batch[blk_idx] = scb2test[j];\
		\
			}\
		\
		}}\


void SymbolicBatchCompressor::compress_symbolic_batch(const astc_codec_image * input_image, const imageblock * blk_batch, symbolic_compressed_block * scb_batch, int cur_batch_size)
{
	cl_int status;
	static_assert((PARTITION_CANDIDATES % 2) == 0, "PARTITION_CANDIDATES should be even number");
	size_t total_finished_blocks = 0;
	batch_size = cur_batch_size;

	OCL_WRITE_BUFFER(blk_buf, sizeof(imageblock) * batch_size, blk_batch);


	memset(blk_stat.host_ptr, BLOCK_STAT_SKIP_ALL, sizeof(uint8_t) * batch_size);

	// compress constant-color blocks
	for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
	{
		const imageblock * blk = &blk_batch[blk_idx];
		symbolic_compressed_block * scb = &scb_batch[blk_idx];

		if (blk->red_min == blk->red_max && blk->green_min == blk->green_max && blk->blue_min == blk->blue_max && blk->alpha_min == blk->alpha_max)
		{
			// detected a constant-color block. Encode as FP16 if using HDR
			scb->error_block = 0;

			if (rgb_force_use_of_hdr)
			{
				scb->block_mode = -1;
				scb->partition_count = 0;
				scb->constant_color[0] = float_to_sf16(blk->orig_data[0], SF_NEARESTEVEN);
				scb->constant_color[1] = float_to_sf16(blk->orig_data[1], SF_NEARESTEVEN);
				scb->constant_color[2] = float_to_sf16(blk->orig_data[2], SF_NEARESTEVEN);
				scb->constant_color[3] = float_to_sf16(blk->orig_data[3], SF_NEARESTEVEN);
			}
			else
			{
				// Encode as UNORM16 if NOT using HDR.
				scb->block_mode = -2;
				scb->partition_count = 0;
				float red = blk->orig_data[0];
				float green = blk->orig_data[1];
				float blue = blk->orig_data[2];
				float alpha = blk->orig_data[3];
				if (red < 0)
					red = 0;
				else if (red > 1)
					red = 1;
				if (green < 0)
					green = 0;
				else if (green > 1)
					green = 1;
				if (blue < 0)
					blue = 0;
				else if (blue > 1)
					blue = 1;
				if (alpha < 0)
					alpha = 0;
				else if (alpha > 1)
					alpha = 1;
				scb->constant_color[0] = (int)floor(red * 65535.0f + 0.5f);
				scb->constant_color[1] = (int)floor(green * 65535.0f + 0.5f);
				scb->constant_color[2] = (int)floor(blue * 65535.0f + 0.5f);
				scb->constant_color[3] = (int)floor(alpha * 65535.0f + 0.5f);
			}

			blk_stat[blk_idx] = BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF;
			total_finished_blocks++;
		}
		else
		{
			error_weight_block *ewb = &ewb_batch[blk_idx];
			error_weight_sum_batch[blk_idx] = prepare_error_weight_block(input_image, xdim, ydim, zdim, &ewp, blk, ewb);
			error_of_best_block[blk_idx] = 1e20f;
		}
	}

	if (total_finished_blocks == batch_size)
		return;

	// compression of average-color blocks disabled for the time being;
	// they produce extremely severe block artifacts.
#if 0
	// first, compress an averaged-color block
	compress_constant_color_block(xdim, ydim, zdim, blk, ewb, scb);

	decompress_symbolic_block(decode_mode, xdim, ydim, zdim, xpos, ypos, zpos, scb, temp);

	float avgblock_errorval = compute_imageblock_difference(xdim, ydim, zdim,
		blk, temp, ewb) * 4.0f;	// bias somewhat against the average-color block.

	if (avgblock_errorval < error_of_best_block)
	{
		error_of_best_block = avgblock_errorval;
		// *scb = scb_candidates[j];
		modesel = 0;
	}

	if ((error_of_best_block / error_weight_sum) < ewp.texel_avg_error_limit)
		goto END_OF_TESTS;
#endif

	ewb_batch.write_to_device();


	float best_errorval_in_mode;
	int i, j;
	uint8_t skip_mode = BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF;


	// next, test mode #0. This mode uses 1 plane of weights and 1 partition.
	compress_symbolic_batch_fixed_partition_1_plane(1, 0, blk_batch, scb_candidates);
		
	for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
	{
		if (blk_stat[blk_idx] & skip_mode)
			continue;

		FIND_BEST_SCB_CANDIDATES(0);
		best_errorvals_in_1pl_1partition_mode[blk_idx] = error_of_best_block[blk_idx];
		if ((error_of_best_block[blk_idx] / error_weight_sum_batch[blk_idx]) < ewp.texel_avg_error_limit)
		{
			blk_stat[blk_idx] |= BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF;
			total_finished_blocks++;
		};
	}

	if (total_finished_blocks == batch_size)
		return;

	//prepare block statistics
	for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
	{
		if (blk_stat[blk_idx] & skip_mode)
			continue;

		const imageblock * blk = &blk_batch[blk_idx];
		const error_weight_block *ewb = &ewb_batch[blk_idx];

		int is_normal_map;
		float lowest_correl;
		
		prepare_block_statistics(xdim, ydim, zdim, blk, ewb, &is_normal_map, &lowest_correl);

		if (is_normal_map)
		{
			blk_stat[blk_idx] |= BLOCK_STAT_NORMAL_MAP;
			if (lowest_correl < 0.99f)
				lowest_correl = 0.99f;
		}
		
		if (lowest_correl > ewp.lowest_correlation_cutoff)
			blk_stat[blk_idx] |= BLOCK_STAT_LOWEST_CORREL_CUTOFF;

		if (blk->alpha_max == blk->alpha_min)
			blk_stat[blk_idx] |= BLOCK_STAT_NO_ALPHA;

		if (blk->grayscale)
			blk_stat[blk_idx] |= BLOCK_STAT_GREYSCALE;

		best_errorvals_in_1pl_2partition_mode[blk_idx] = 1e30f;
	}

	// next, test the four possible 1-partition, 2-planes modes
	skip_mode = BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF | BLOCK_STAT_LOWEST_CORREL_CUTOFF | BLOCK_STAT_GREYSCALE;
	for (i = 0; i < 4; i++)
	{
		if (i == 3)
			skip_mode = BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF | BLOCK_STAT_LOWEST_CORREL_CUTOFF | BLOCK_STAT_NO_ALPHA;

		compress_symbolic_batch_fixed_partition_2_planes(1, 0, i, blk_batch, scb_candidates, skip_mode);

		for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
		{
			if (blk_stat[blk_idx] & skip_mode)
				continue;

			FIND_BEST_SCB_CANDIDATES(i + 1);
			if ((error_of_best_block[blk_idx] / error_weight_sum_batch[blk_idx]) < ewp.texel_avg_error_limit)
			{
				blk_stat[blk_idx] |= BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF;
				total_finished_blocks++;
			};
		}

		if (total_finished_blocks == batch_size)
			return;
	}



	// find best blocks for 2, 3 and 4 partitions
	for (int partition_count = 2; partition_count <= 4; partition_count++)
	{
		//find_best_partitionings_batch_ocl(partition_count, blk_batch);
		find_best_partitionings_batch(partition_count, blk_batch);

		skip_mode = BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF;
		for (i = 0; i < 2; i++)
		{
			compress_symbolic_batch_fixed_partition_1_plane(partition_count, i, blk_batch, scb_candidates);
			
			for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
			{
				if (blk_stat[blk_idx] & skip_mode)
					continue;

				FIND_BEST_SCB_CANDIDATES(4 * (partition_count - 2) + 5 + i);
				best_errorvals_in_1pl_2partition_mode[blk_idx] = MIN(best_errorval_in_mode, best_errorvals_in_1pl_2partition_mode[blk_idx]);
				if ((error_of_best_block[blk_idx] / error_weight_sum_batch[blk_idx]) < ewp.texel_avg_error_limit)
				{
					blk_stat[blk_idx] |= BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF;
					total_finished_blocks++;
				};
			}
		}

		// cut compression work off if haven't gain much from 2 partition mode
		if (partition_count == 2)
		for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
		{
			if (blk_stat[blk_idx] & (BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF | BLOCK_STAT_NORMAL_MAP))
				continue;

			if (best_errorvals_in_1pl_2partition_mode[blk_idx] > (best_errorvals_in_1pl_1partition_mode[blk_idx] * ewp.partition_1_to_2_limit))
			{
				blk_stat[blk_idx] |= BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF;
				total_finished_blocks++;
			}
		}

		// don't bother to check 4 partitions for dual plane of weightss, ever.
		if (partition_count == 4)
			return;

		if (total_finished_blocks == batch_size)
			return;

		skip_mode = BLOCK_STAT_LOWEST_CORREL_CUTOFF | BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF;
		for (i = 0; i < 2; i++)
		{
			compress_symbolic_batch_fixed_partition_2_planes(partition_count, i, 0, blk_batch, scb_candidates, skip_mode);
			
			for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
			{
				if (blk_stat[blk_idx] & skip_mode)
					continue;

				FIND_BEST_SCB_CANDIDATES(4 * (partition_count - 2) + 5 + 2 + i);
				if ((error_of_best_block[blk_idx] / error_weight_sum_batch[blk_idx]) < ewp.texel_avg_error_limit)
				{
					blk_stat[blk_idx] |= BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF;
					total_finished_blocks++;
				};
			}
		}
	}
}

SymbolicBatchCompressor::~SymbolicBatchCompressor()
{
	cl_int status;

	delete[] tmpplanes.per_scb.u8_qdq_weights;
	delete[] tmpplanes.per_scb.decimation_mode;
	delete[] tmpplanes.per_scb.color_quantization_level_mod;
	delete[] tmpplanes.per_scb.color_quantization_level;
	delete[] tmpplanes.per_scb.weight_mode;
	delete[] tmpplanes.per_scb.partition_format_specifiers;
	delete[] tmpplanes.scb_stat;
	delete[] tmpplanes.weight_high_value2;
	delete[] tmpplanes.weight_high_value1;
	delete[] tmpplanes.weight_low_value2;
	delete[] tmpplanes.weight_low_value1;
	delete[] tmpplanes.u8_quantized_decimated_quantized_weights;
	delete[] tmpplanes.flt_quantized_decimated_quantized_weights;
	delete[] tmpplanes.decimated_weights;
	delete[] tmpplanes.decimated_quantized_weights;
	delete[] tmpplanes.per_scb.ep;
	delete[] tmpplanes.eix2;
	delete[] tmpplanes.eix1;
	delete[] tmpplanes.ei2;
	delete[] tmpplanes.ei1;

	delete[] error_of_best_block;
	delete[] best_errorvals_in_1pl_2partition_mode;
	delete[] best_errorvals_in_1pl_1partition_mode;
	delete[] error_weight_sum_batch;
	delete[] scb_candidates;

	OCL_RELEASE_OBJECT(Kernel, fbp.find_best_partitionings);
	for (size_t pcount = 4; pcount >= 2; pcount--)
	{
		OCL_RELEASE_OBJECT(MemObject, fbp.ptab[pcount]);
	}
	OCL_RELEASE_OBJECT(MemObject, blk_buf);

	OCL_RELEASE_OBJECT(CommandQueue, opencl_queue);
}

void SymbolicBatchCompressor::allocate_buffers(int max_blocks)
{
	ewb_batch.create_buffer(CL_MEM_READ_ONLY, max_blocks, opencl_queue);
	scb_candidates = new symbolic_compressed_block[SCB_CANDIDATES * max_blocks];
	blk_stat.create_buffer(CL_MEM_READ_ONLY, max_blocks, opencl_queue);
	error_weight_sum_batch = new float[max_blocks];
	best_errorvals_in_1pl_1partition_mode = new float[max_blocks];
	best_errorvals_in_1pl_2partition_mode = new float[max_blocks];
	error_of_best_block = new float[max_blocks];

	partition_indices_1plane_batch.create_buffer(CL_MEM_WRITE_ONLY, PARTITION_CANDIDATES * max_blocks, opencl_queue);
	partition_indices_2planes_batch.create_buffer(CL_MEM_WRITE_ONLY, PARTITION_CANDIDATES * max_blocks, opencl_queue);

	tmpplanes.ei1 = new endpoints_and_weights[max_blocks];
	tmpplanes.ei2 = new endpoints_and_weights[max_blocks];
	tmpplanes.eix1 = new endpoints_and_weights[MAX_DECIMATION_MODES * max_blocks];
	tmpplanes.eix2 = new endpoints_and_weights[MAX_DECIMATION_MODES * max_blocks];
	tmpplanes.per_scb.ep = new endpoints[SCB_CANDIDATES * max_blocks];
	tmpplanes.decimated_quantized_weights = new float[MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK * max_blocks];
	tmpplanes.decimated_weights = new float[MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK * max_blocks];
	tmpplanes.flt_quantized_decimated_quantized_weights = new float[MAX_SORTED_WEIGHT_MODES * MAX_WEIGHTS_PER_BLOCK * max_blocks];
	tmpplanes.u8_quantized_decimated_quantized_weights = new uint8_t[MAX_SORTED_WEIGHT_MODES * MAX_WEIGHTS_PER_BLOCK * max_blocks];
	tmpplanes.weight_low_value1 = new float[MAX_SORTED_WEIGHT_MODES * max_blocks];
	tmpplanes.weight_low_value2 = new float[MAX_SORTED_WEIGHT_MODES * max_blocks];
	tmpplanes.weight_high_value1 = new float[MAX_SORTED_WEIGHT_MODES * max_blocks];
	tmpplanes.weight_high_value2 = new float[MAX_SORTED_WEIGHT_MODES * max_blocks];
	tmpplanes.scb_stat = new uint8_t[max_blocks];
	tmpplanes.per_scb.partition_format_specifiers = new int[SCB_CANDIDATES * MAX_PARTITIONS * max_blocks];
	tmpplanes.per_scb.weight_mode = new int[SCB_CANDIDATES * max_blocks];
	tmpplanes.per_scb.color_quantization_level = new int[SCB_CANDIDATES * max_blocks];
	tmpplanes.per_scb.color_quantization_level_mod = new int[SCB_CANDIDATES * max_blocks];
	tmpplanes.per_scb.decimation_mode = new int[SCB_CANDIDATES * max_blocks];
	tmpplanes.per_scb.u8_qdq_weights = new uint8_t[MAX_WEIGHTS_PER_BLOCK * SCB_CANDIDATES * max_blocks];

	fbp.partition_sequence.create_buffer(CL_MEM_READ_ONLY, max_blocks * PARTITION_COUNT, opencl_queue);
}


// function for compressing a block symbolically, given
// that we have already decided on a partition
void SymbolicBatchCompressor::compress_symbolic_batch_fixed_partition_1_plane(int partition_count, int partition_offset, const imageblock * blk_batch, symbolic_compressed_block * scb_candidates)
{
	static const int free_bits_for_partition_count_1plane[5] = { 0, 115 - 4, 111 - 4 - PARTITION_BITS, 108 - 4 - PARTITION_BITS, 105 - 4 - PARTITION_BITS };
	// maximum bits available for weights = (ei_bits - color_bits)
	// ei_bits - number of bits available for weights and color data if all color endpoints are of the same type
	// color_bits - minimum number of bits needed to store color data (QUANT_5, FMT_LUMINANCE)
	static const int max_weight_bits_for_partition_count_1plane[5] = { 0, 111 - 5, 99 - 10, 99 - 14, 99 - 19 };
	float mode_cutoff = ewp.block_mode_cutoff;

	const block_size_descriptor_sorted *sorted_bsd = get_sorted_block_size_descriptor(xdim, ydim, zdim, 0);
	const decimation_table *const *ixtab3 = sorted_bsd->decimation_tables;
	const partition_info *ptab = get_partition_table(xdim, ydim, zdim, partition_count);


	for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
	{
		if (blk_stat[blk_idx] & BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF)
			continue;

		const imageblock * blk = &blk_batch[blk_idx];
		const error_weight_block *ewb = &ewb_batch[blk_idx];
		endpoints_and_weights *ei = &tmpplanes.ei1[blk_idx];
		endpoints_and_weights *eix = &tmpplanes.eix1[MAX_DECIMATION_MODES * blk_idx];
		float *decimated_quantized_weights = &tmpplanes.decimated_quantized_weights[MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];
		float *decimated_weights = &tmpplanes.decimated_weights[MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];

		int partition_index = 0;
		if (partition_count > 1)
		{
			partition_index = partition_indices_1plane_batch[blk_idx * PARTITION_CANDIDATES + partition_offset];
		}
		const partition_info *pi = &ptab[partition_index];

		// first, compute ideal weights and endpoint colors, under thre assumption that
		// there is no quantization or decimation going on.
		compute_endpoints_and_ideal_weights_1_plane(xdim, ydim, zdim, pi, blk, ewb, ei);


		// next, compute ideal weights and endpoint colors for every decimation.

		// for each decimation mode, compute an ideal set of weights
		// (that is, weights computed with the assumption that they are not quantized)
		for (int i = 0; i < ewp.decimation_mode_limit_1plane; i++)
		{
			eix[i] = *ei;
			compute_ideal_weights_for_decimation_table(&(eix[i]), ixtab3[i], decimated_quantized_weights + i * MAX_WEIGHTS_PER_BLOCK, decimated_weights + i * MAX_WEIGHTS_PER_BLOCK);
		}
	}


	// for each mode, use the angular method to compute a shift.
	for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
	{
		if (blk_stat[blk_idx] & BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF)
			continue;

		const float *decimated_quantized_weights = &tmpplanes.decimated_quantized_weights[MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];
		const float *decimated_weights = &tmpplanes.decimated_weights[MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];
		float *weight_low_value = &tmpplanes.weight_low_value1[MAX_SORTED_WEIGHT_MODES * blk_idx];
		float *weight_high_value = &tmpplanes.weight_high_value1[MAX_SORTED_WEIGHT_MODES * blk_idx];

		compute_angular_endpoints_1plane(&ewp, sorted_bsd, decimated_quantized_weights, decimated_weights, weight_low_value, weight_high_value);
	}

	memset(tmpplanes.scb_stat, 0, sizeof(uint8_t) * batch_size);
	for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
	{
		if (blk_stat[blk_idx] & BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF)
		{
			tmpplanes.scb_stat[blk_idx] = 0x0F;
			continue;
		}

		const imageblock * blk = &blk_batch[blk_idx];
		const error_weight_block *ewb = &ewb_batch[blk_idx];
		symbolic_compressed_block *scb = &scb_candidates[SCB_CANDIDATES * blk_idx];
		const endpoints_and_weights *ei = &tmpplanes.ei1[blk_idx];
		endpoints_and_weights *eix = &tmpplanes.eix1[MAX_DECIMATION_MODES * blk_idx];
		const float *decimated_quantized_weights = &tmpplanes.decimated_quantized_weights[MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];
		float *flt_quantized_decimated_quantized_weights = &tmpplanes.flt_quantized_decimated_quantized_weights[MAX_SORTED_WEIGHT_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];
		uint8_t *u8_quantized_decimated_quantized_weights = &tmpplanes.u8_quantized_decimated_quantized_weights[MAX_SORTED_WEIGHT_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];
		const float *weight_low_value = &tmpplanes.weight_low_value1[MAX_SORTED_WEIGHT_MODES * blk_idx];
		float *weight_high_value = &tmpplanes.weight_high_value1[MAX_SORTED_WEIGHT_MODES * blk_idx];

		int partition_index = 0;
		if (partition_count > 1)
		{
			partition_index = partition_indices_1plane_batch[blk_idx * PARTITION_CANDIDATES + partition_offset];
		}
		const partition_info *pi = &ptab[partition_index];


		// compute maximum colors for the endpoints and ideal weights.
		// for each endpoint-and-ideal-weight pair, compute the smallest weight value
		// that will result in a color value greater than 1.
		float4 min_ep = float4(10, 10, 10, 10);
		for (int i = 0; i < partition_count; i++)
		{
			float4 ep = (float4(1, 1, 1, 1) - ei->ep.endpt0[i]) / (ei->ep.endpt1[i] - ei->ep.endpt0[i]);
			if (ep.x > 0.5f && ep.x < min_ep.x)
				min_ep.x = ep.x;
			if (ep.y > 0.5f && ep.y < min_ep.y)
				min_ep.y = ep.y;
			if (ep.z > 0.5f && ep.z < min_ep.z)
				min_ep.z = ep.z;
			if (ep.w > 0.5f && ep.w < min_ep.w)
				min_ep.w = ep.w;
		}
		float min_wt_cutoff = MIN(MIN(min_ep.x, min_ep.y), MIN(min_ep.z, min_ep.w));

		// for each mode (which specifies a decimation and a quantization):
		// * compute number of bits needed for the quantized weights.
		// * generate an optimized set of quantized weights.
		// * compute quantization errors for the mode.

		int qwt_bitcounts[MAX_SORTED_WEIGHT_MODES];
		float qwt_errors[MAX_SORTED_WEIGHT_MODES];
		for (int i = 0; i < ewp.weight_mode_limit_1plane; i++)
		{
			if (weight_high_value[i] > 1.02f * min_wt_cutoff)
				weight_high_value[i] = 1.0f;

			int decimation_mode = sorted_bsd->block_modes[i].sorted_decimation_mode;


			// compute weight bitcount for the mode
			int bits_used_by_weights = compute_ise_bitcount(ixtab3[decimation_mode]->num_weights,
				(quantization_method)sorted_bsd->block_modes[i].quantization_mode);
			int bitcount = free_bits_for_partition_count_1plane[partition_count] - bits_used_by_weights;
			if (bits_used_by_weights > max_weight_bits_for_partition_count_1plane[partition_count])
			{
				qwt_errors[i] = 1e38f;
				continue;
			}
			qwt_bitcounts[i] = bitcount;


			// then, generate the optimized set of weights for the weight mode.
			compute_ideal_quantized_weights_for_decimation_table(&(eix[decimation_mode]),
				ixtab3[decimation_mode],
				weight_low_value[i], weight_high_value[i],
				decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * decimation_mode,
				flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * i,
				u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * i,
				sorted_bsd->block_modes[i].quantization_mode);

			// then, compute weight-errors for the weight mode.
			qwt_errors[i] = compute_error_of_weight_set(&(eix[decimation_mode]), ixtab3[decimation_mode], flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * i);
		}

		// for each weighting mode, determine the optimal combination of color endpoint encodings
		// and weight encodings; return results for the 4 best-looking modes.

		int *partition_format_specifiers = &tmpplanes.per_scb.partition_format_specifiers[MAX_PARTITIONS * SCB_CANDIDATES * blk_idx];
		int *weight_mode = &tmpplanes.per_scb.weight_mode[SCB_CANDIDATES * blk_idx];
		int *color_quantization_level = &tmpplanes.per_scb.color_quantization_level[SCB_CANDIDATES * blk_idx];
		int *color_quantization_level_mod = &tmpplanes.per_scb.color_quantization_level_mod[SCB_CANDIDATES * blk_idx];
		determine_optimal_set_of_endpoint_formats_to_use(xdim, ydim, zdim, pi, blk, ewb, &(ei->ep), -1,	// used to flag that we are in single-weight mode
			qwt_bitcounts, qwt_errors, ewp.weight_mode_limit_1plane, partition_format_specifiers, weight_mode, color_quantization_level, color_quantization_level_mod);

		for (int i = 0; i < 4; i++)
		{
			if (weight_mode[i] < 0)
			{
				scb[i].error_block = 1;
				tmpplanes.scb_stat[blk_idx] |= 1 << i;
				continue;
			}

			int scb_idx = blk_idx * SCB_CANDIDATES + i;
			int decimation_mode = sorted_bsd->block_modes[weight_mode[i]].sorted_decimation_mode;
			tmpplanes.per_scb.decimation_mode[scb_idx] = decimation_mode;
			tmpplanes.per_scb.ep[scb_idx] = eix[decimation_mode].ep;

			const uint8_t *u8_weight_src = u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * weight_mode[i];
			memcpy(&tmpplanes.per_scb.u8_qdq_weights[MAX_WEIGHTS_PER_BLOCK * scb_idx], u8_weight_src, sizeof(uint8_t) * MAX_WEIGHTS_PER_BLOCK);
		}
	}

	
	// then iterate over the 4 believed-to-be-best modes to find out which one is
	// actually best.
	for (int l = 0; l < ewp.max_refinement_iters; l++)
	{
		for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
		{
			if (blk_stat[blk_idx] & BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF)
				continue;

			const imageblock * blk = &blk_batch[blk_idx];
			const error_weight_block *ewb = &ewb_batch[blk_idx];
			symbolic_compressed_block *scb = &scb_candidates[SCB_CANDIDATES * blk_idx];

			int partition_index = 0;
			if (partition_count > 1)
			{
				partition_index = partition_indices_1plane_batch[blk_idx * PARTITION_CANDIDATES + partition_offset];
			}
			const partition_info *pi = &ptab[partition_index];

			for (int i = 0; i < 4; i++)
			{
				if (tmpplanes.scb_stat[blk_idx] & (0x11 << i))
					continue;

				int scb_idx = SCB_CANDIDATES * blk_idx + i;
				
				const int weight_mode = tmpplanes.per_scb.weight_mode[scb_idx];
				const int color_quantization_level = tmpplanes.per_scb.color_quantization_level[scb_idx];
				const int color_quantization_level_mod = tmpplanes.per_scb.color_quantization_level_mod[scb_idx];
				const int decimation_mode = tmpplanes.per_scb.decimation_mode[scb_idx];
				const int *partition_format_specifiers = &tmpplanes.per_scb.partition_format_specifiers[MAX_PARTITIONS * scb_idx];
				const int weight_quantization_mode = sorted_bsd->block_modes[weight_mode].quantization_mode;
				const decimation_table *it = ixtab3[decimation_mode];

				const uint8_t *u8_weight_src = &tmpplanes.per_scb.u8_qdq_weights[MAX_WEIGHTS_PER_BLOCK * scb_idx];

				// recompute the ideal color endpoints before storing them.
				endpoints *ep = &tmpplanes.per_scb.ep[scb_idx];
				float4 rgbs_colors[4];
				float4 rgbo_colors[4];
				float2 lum_intervals[4];


				recompute_ideal_colors(xdim, ydim, zdim, weight_quantization_mode, ep, rgbs_colors, rgbo_colors, lum_intervals, u8_weight_src, NULL, -1, pi, it, blk, ewb);

				// quantize the chosen color

				// store the colors for the block
				for (int j = 0; j < partition_count; j++)
				{
					scb[i].color_formats[j] = pack_color_endpoints(decode_mode,
						ep->endpt0[j],
						ep->endpt1[j],
						rgbs_colors[j], rgbo_colors[j], lum_intervals[j], partition_format_specifiers[j], scb[i].color_values[j], color_quantization_level);
				}


				// if all the color endpoint modes are the same, we get a few more
				// bits to store colors; let's see if we can take advantage of this:
				// requantize all the colors and see if the endpoint modes remain the same;
				// if they do, then exploit it.
				scb[i].color_formats_matched = 0;

				if ((partition_count >= 2 && scb[i].color_formats[0] == scb[i].color_formats[1]
					&& color_quantization_level != color_quantization_level_mod)
					&& (partition_count == 2 || (scb[i].color_formats[0] == scb[i].color_formats[2] && (partition_count == 3 || (scb[i].color_formats[0] == scb[i].color_formats[3])))))
				{
					int colorvals[4][12];
					int color_formats_mod[4];
					for (int j = 0; j < partition_count; j++)
					{
						color_formats_mod[j] = pack_color_endpoints(decode_mode,
							ep->endpt0[j],
							ep->endpt1[j],
							rgbs_colors[j], rgbo_colors[j], lum_intervals[j], partition_format_specifiers[j], colorvals[j], color_quantization_level_mod);
					}
					if (color_formats_mod[0] == color_formats_mod[1]
						&& (partition_count == 2 || (color_formats_mod[0] == color_formats_mod[2] && (partition_count == 3 || (color_formats_mod[0] == color_formats_mod[3])))))
					{
						scb[i].color_formats_matched = 1;
						for (int j = 0; j < 4; j++)
							for (int k = 0; k < 12; k++)
								scb[i].color_values[j][k] = colorvals[j][k];
						for (int j = 0; j < 4; j++)
							scb[i].color_formats[j] = color_formats_mod[j];
					}
				}


				// store header fields
				scb[i].partition_count = partition_count;
				scb[i].partition_index = partition_index;
				scb[i].color_quantization_level = scb[i].color_formats_matched ? color_quantization_level_mod : color_quantization_level;
				scb[i].block_mode = sorted_bsd->block_modes[weight_mode].block_mode;
				scb[i].error_block = 0;

				if (scb[i].color_quantization_level < 4)
				{
					scb[i].error_block = 1;	// should never happen, but cannot prove it impossible.
				}

			}
		}

		for (int scb_idx = 0; scb_idx < batch_size * SCB_CANDIDATES; scb_idx++)
		{
			int blk_idx = scb_idx >> 2;
			int scb_off = scb_idx & 0x03;
			if (tmpplanes.scb_stat[blk_idx] & (0x11 << scb_off))
				continue;

			const imageblock * blk = &blk_batch[blk_idx];
			const error_weight_block *ewb = &ewb_batch[blk_idx];

			uint8_t *u8_weight_src = &tmpplanes.per_scb.u8_qdq_weights[MAX_WEIGHTS_PER_BLOCK * scb_idx];

			int adjustments = realign_weights(decode_mode,
				xdim, ydim, zdim,
				blk, ewb, &scb_candidates[scb_idx],
				u8_weight_src,
				NULL);

			if (adjustments == 0)
				tmpplanes.scb_stat[blk_idx] |= 0x10 << scb_off;
		}
	}

	for (int scb_idx = 0; scb_idx < batch_size * SCB_CANDIDATES; scb_idx++)
	{
		int blk_idx = scb_idx >> 2;
		int scb_off = scb_idx & 0x03;
		if (tmpplanes.scb_stat[blk_idx] & (1 << scb_off))
			continue;
		
		int decimation_mode = tmpplanes.per_scb.decimation_mode[scb_idx];
		const decimation_table *it = ixtab3[decimation_mode];

		const uint8_t *u8_weight_src = &tmpplanes.per_scb.u8_qdq_weights[MAX_WEIGHTS_PER_BLOCK * scb_idx];
		int weights_to_copy = it->num_weights;

		for (int j = 0; j < weights_to_copy; j++)
			scb_candidates[scb_idx].plane1_weights[j] = u8_weight_src[j];
	}
}


void SymbolicBatchCompressor::compress_symbolic_batch_fixed_partition_2_planes(int partition_count, int partition_offset, int separate_component, const imageblock * blk_batch, symbolic_compressed_block * scb_candidates, uint8_t skip_mode)
{
	static const int free_bits_for_partition_count_2planes[5] = { 0, 113 - 4, 109 - 4 - PARTITION_BITS, 106 - 4 - PARTITION_BITS, 103 - 4 - PARTITION_BITS };
	// maximum bits available for weights = (ei_bits - color_bits)
	// ei_bits - number of bits available for weights and color data if all color endpoints are of the same type
	// color_bits - minimum number of bits needed to store color data (QUANT_5, FMT_LUMINANCE)
	static const int max_weight_bits_for_partition_count_2planes[5] = { 0, 109 - 5, 97 - 10, 97 - 14, 97 - 19 };
	float mode_cutoff = ewp.block_mode_cutoff;

	const block_size_descriptor_sorted *sorted_bsd = get_sorted_block_size_descriptor(xdim, ydim, zdim, 1);
	const decimation_table *const *ixtab3 = sorted_bsd->decimation_tables;
	const partition_info *ptab = get_partition_table(xdim, ydim, zdim, partition_count);


	for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
	{
		if (blk_stat[blk_idx] & skip_mode)
			continue;

		const imageblock * blk = &blk_batch[blk_idx];
		const error_weight_block *ewb = &ewb_batch[blk_idx];
		endpoints_and_weights *ei1 = &tmpplanes.ei1[blk_idx];
		endpoints_and_weights *ei2 = &tmpplanes.ei2[blk_idx];
		endpoints_and_weights *eix1 = &tmpplanes.eix1[MAX_DECIMATION_MODES * blk_idx];
		endpoints_and_weights *eix2 = &tmpplanes.eix2[MAX_DECIMATION_MODES * blk_idx];
		float *decimated_quantized_weights = &tmpplanes.decimated_quantized_weights[MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];
		float *decimated_weights = &tmpplanes.decimated_weights[MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];

		int partition_index = 0;
		if (partition_count > 1)
		{
			partition_index = partition_indices_2planes_batch[blk_idx * PARTITION_CANDIDATES + partition_offset];
			separate_component = partition_index >> PARTITION_BITS;
			partition_index = partition_index & (PARTITION_COUNT - 1);
		}
		const partition_info *pi = &ptab[partition_index];

		// first, compute ideal weights and endpoint colors
		compute_endpoints_and_ideal_weights_2_planes(xdim, ydim, zdim, pi, blk, ewb, separate_component, ei1, ei2);


		// next, compute ideal weights and endpoint colors for every decimation.

		// for each decimation mode, compute an ideal set of weights
		for (int i = 0; i < ewp.decimation_mode_limit_2planes; i++)
		{
			eix1[i] = *ei1;
			eix2[i] = *ei2;
			compute_ideal_weights_for_decimation_table(&(eix1[i]), ixtab3[i], decimated_quantized_weights + (2 * i) * (MAX_WEIGHTS_PER_BLOCK / 2), decimated_weights + (2 * i) * (MAX_WEIGHTS_PER_BLOCK / 2));
			compute_ideal_weights_for_decimation_table(&(eix2[i]), ixtab3[i], decimated_quantized_weights + (2 * i + 1) * (MAX_WEIGHTS_PER_BLOCK / 2), decimated_weights + (2 * i + 1) * (MAX_WEIGHTS_PER_BLOCK / 2));
		}
	}


	// for each mode, use the angular method to compute a shift.
	for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
	{
		if (blk_stat[blk_idx] & skip_mode)
			continue;

		const float *decimated_quantized_weights = &tmpplanes.decimated_quantized_weights[MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];
		const float *decimated_weights = &tmpplanes.decimated_weights[MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];
		float *weight_low_value1 = &tmpplanes.weight_low_value1[MAX_SORTED_WEIGHT_MODES * blk_idx];
		float *weight_low_value2 = &tmpplanes.weight_low_value2[MAX_SORTED_WEIGHT_MODES * blk_idx];
		float *weight_high_value1 = &tmpplanes.weight_high_value1[MAX_SORTED_WEIGHT_MODES * blk_idx];
		float *weight_high_value2 = &tmpplanes.weight_high_value2[MAX_SORTED_WEIGHT_MODES * blk_idx];

		compute_angular_endpoints_2planes(&ewp, sorted_bsd, decimated_quantized_weights, decimated_weights, weight_low_value1, weight_high_value1, weight_low_value2, weight_high_value2);
	}


	memset(tmpplanes.scb_stat, 0, sizeof(uint8_t) * batch_size);
	for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
	{
		if (blk_stat[blk_idx] & skip_mode)
		{
			tmpplanes.scb_stat[blk_idx] = 0x0F;
			continue;
		}

		const imageblock * blk = &blk_batch[blk_idx];
		const error_weight_block *ewb = &ewb_batch[blk_idx];
		symbolic_compressed_block *scb = &scb_candidates[SCB_CANDIDATES * blk_idx];
		const endpoints_and_weights *ei1 = &tmpplanes.ei1[blk_idx];
		const endpoints_and_weights *ei2 = &tmpplanes.ei2[blk_idx];
		endpoints_and_weights *eix1 = &tmpplanes.eix1[MAX_DECIMATION_MODES * blk_idx];
		endpoints_and_weights *eix2 = &tmpplanes.eix2[MAX_DECIMATION_MODES * blk_idx];
		const float *decimated_quantized_weights = &tmpplanes.decimated_quantized_weights[MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];
		float *flt_quantized_decimated_quantized_weights = &tmpplanes.flt_quantized_decimated_quantized_weights[MAX_SORTED_WEIGHT_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];
		uint8_t *u8_quantized_decimated_quantized_weights = &tmpplanes.u8_quantized_decimated_quantized_weights[MAX_SORTED_WEIGHT_MODES * MAX_WEIGHTS_PER_BLOCK * blk_idx];
		const float *weight_low_value1 = &tmpplanes.weight_low_value1[MAX_SORTED_WEIGHT_MODES * blk_idx];
		const float *weight_low_value2 = &tmpplanes.weight_low_value2[MAX_SORTED_WEIGHT_MODES * blk_idx];
		float *weight_high_value1 = &tmpplanes.weight_high_value1[MAX_SORTED_WEIGHT_MODES * blk_idx];
		float *weight_high_value2 = &tmpplanes.weight_high_value2[MAX_SORTED_WEIGHT_MODES * blk_idx];

		int partition_index = 0;
		if (partition_count > 1)
		{
			partition_index = partition_indices_2planes_batch[blk_idx * PARTITION_CANDIDATES + partition_offset];
			separate_component = partition_index >> PARTITION_BITS;
			partition_index = partition_index & (PARTITION_COUNT - 1);
		}
		const partition_info *pi = &ptab[partition_index];

		// compute maximum colors for the endpoints and ideal weights.
		// for each endpoint-and-ideal-weight pair, compute the smallest weight value
		// that will result in a color value greater than 1.
		float4 min_ep1 = float4(10, 10, 10, 10);
		float4 min_ep2 = float4(10, 10, 10, 10);
		for (int i = 0; i < partition_count; i++)
		{
			float4 ep1 = (float4(1, 1, 1, 1) - ei1->ep.endpt0[i]) / (ei1->ep.endpt1[i] - ei1->ep.endpt0[i]);
			if (ep1.x > 0.5f && ep1.x < min_ep1.x)
				min_ep1.x = ep1.x;
			if (ep1.y > 0.5f && ep1.y < min_ep1.y)
				min_ep1.y = ep1.y;
			if (ep1.z > 0.5f && ep1.z < min_ep1.z)
				min_ep1.z = ep1.z;
			if (ep1.w > 0.5f && ep1.w < min_ep1.w)
				min_ep1.w = ep1.w;
			float4 ep2 = (float4(1, 1, 1, 1) - ei2->ep.endpt0[i]) / (ei2->ep.endpt1[i] - ei2->ep.endpt0[i]);
			if (ep2.x > 0.5f && ep2.x < min_ep2.x)
				min_ep2.x = ep2.x;
			if (ep2.y > 0.5f && ep2.y < min_ep2.y)
				min_ep2.y = ep2.y;
			if (ep2.z > 0.5f && ep2.z < min_ep2.z)
				min_ep2.z = ep2.z;
			if (ep2.w > 0.5f && ep2.w < min_ep2.w)
				min_ep2.w = ep2.w;
		}

		float min_wt_cutoff1, min_wt_cutoff2;
		switch (separate_component)
		{
		case 0:
			min_wt_cutoff2 = min_ep2.x;
			min_ep1.x = 1e30f;
			break;
		case 1:
			min_wt_cutoff2 = min_ep2.y;
			min_ep1.y = 1e30f;
			break;
		case 2:
			min_wt_cutoff2 = min_ep2.z;
			min_ep1.z = 1e30f;
			break;
		case 3:
			min_wt_cutoff2 = min_ep2.w;
			min_ep1.w = 1e30f;
			break;
		default:
			min_wt_cutoff2 = 1e30f;
		}

		min_wt_cutoff1 = MIN(MIN(min_ep1.x, min_ep1.y), MIN(min_ep1.z, min_ep1.w));

		// for each mode (which specifies a decimation and a quantization):
		// * compute number of bits needed for the quantized weights.
		// * generate an optimized set of quantized weights.
		// * compute quantization errors for each mode

		int qwt_bitcounts[MAX_SORTED_WEIGHT_MODES];
		float qwt_errors[MAX_SORTED_WEIGHT_MODES];
		for (int i = 0; i < ewp.weight_mode_limit_2planes; i++)
		{
			int decimation_mode = sorted_bsd->block_modes[i].sorted_decimation_mode;

			if (weight_high_value1[i] > 1.02f * min_wt_cutoff1)
				weight_high_value1[i] = 1.0f;
			if (weight_high_value2[i] > 1.02f * min_wt_cutoff2)
				weight_high_value2[i] = 1.0f;

			// compute weight bitcount for the mode
			int bits_used_by_weights = compute_ise_bitcount(2 * ixtab3[decimation_mode]->num_weights,
				(quantization_method)sorted_bsd->block_modes[i].quantization_mode);
			int bitcount = free_bits_for_partition_count_2planes[partition_count] - bits_used_by_weights;
			if (bits_used_by_weights > max_weight_bits_for_partition_count_2planes[partition_count])
			{
				qwt_errors[i] = 1e38f;
				continue;
			}
			qwt_bitcounts[i] = bitcount;


			// then, generate the optimized set of weights for the mode.
			compute_ideal_quantized_weights_for_decimation_table(&(eix1[decimation_mode]),
				ixtab3[decimation_mode],
				weight_low_value1[i],
				weight_high_value1[i],
				decimated_quantized_weights + (MAX_WEIGHTS_PER_BLOCK / 2) * (2 * decimation_mode),
				flt_quantized_decimated_quantized_weights + (MAX_WEIGHTS_PER_BLOCK / 2) * (2 * i),
				u8_quantized_decimated_quantized_weights + (MAX_WEIGHTS_PER_BLOCK / 2) * (2 * i), sorted_bsd->block_modes[i].quantization_mode);
			compute_ideal_quantized_weights_for_decimation_table(&(eix2[decimation_mode]),
				ixtab3[decimation_mode],
				weight_low_value2[i],
				weight_high_value2[i],
				decimated_quantized_weights + (MAX_WEIGHTS_PER_BLOCK / 2) * (2 * decimation_mode + 1),
				flt_quantized_decimated_quantized_weights + (MAX_WEIGHTS_PER_BLOCK / 2) * (2 * i + 1),
				u8_quantized_decimated_quantized_weights + (MAX_WEIGHTS_PER_BLOCK / 2) * (2 * i + 1), sorted_bsd->block_modes[i].quantization_mode);


			// then, compute quantization errors for the block mode.
			qwt_errors[i] =
				compute_error_of_weight_set(&(eix1[decimation_mode]),
					ixtab3[decimation_mode],
					flt_quantized_decimated_quantized_weights + (MAX_WEIGHTS_PER_BLOCK / 2) * (2 * i))
				+ compute_error_of_weight_set(&(eix2[decimation_mode]), ixtab3[decimation_mode], flt_quantized_decimated_quantized_weights + (MAX_WEIGHTS_PER_BLOCK / 2) * (2 * i + 1));
		}


		// decide the optimal combination of color endpoint encodings and weight encoodings.
		int *partition_format_specifiers = &tmpplanes.per_scb.partition_format_specifiers[MAX_PARTITIONS * SCB_CANDIDATES * blk_idx];
		int *weight_mode = &tmpplanes.per_scb.weight_mode[SCB_CANDIDATES * blk_idx];
		int *color_quantization_level = &tmpplanes.per_scb.color_quantization_level[SCB_CANDIDATES * blk_idx];
		int *color_quantization_level_mod = &tmpplanes.per_scb.color_quantization_level_mod[SCB_CANDIDATES * blk_idx];

		endpoints epm;
		merge_endpoints(&(ei1->ep), &(ei2->ep), separate_component, &epm);

		determine_optimal_set_of_endpoint_formats_to_use(xdim, ydim, zdim,
			pi,
			blk,
			ewb,
			&epm, separate_component, qwt_bitcounts, qwt_errors, ewp.weight_mode_limit_2planes, partition_format_specifiers, weight_mode, color_quantization_level, color_quantization_level_mod);

		for (int i = 0; i < 4; i++)
		{
			if (weight_mode[i] < 0)
			{
				scb[i].error_block = 1;
				tmpplanes.scb_stat[blk_idx] |= 1 << i;
				continue;
			}

			int scb_idx = blk_idx * SCB_CANDIDATES + i;
			int decimation_mode = sorted_bsd->block_modes[weight_mode[i]].sorted_decimation_mode;
			tmpplanes.per_scb.decimation_mode[scb_idx] = decimation_mode;
			merge_endpoints(&(eix1[decimation_mode].ep), &(eix2[decimation_mode].ep), separate_component, &tmpplanes.per_scb.ep[scb_idx]);

			// copy weights of both planes
			const uint8_t *u8_weight_src = u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * weight_mode[i];
			memcpy(&tmpplanes.per_scb.u8_qdq_weights[MAX_WEIGHTS_PER_BLOCK * scb_idx], u8_weight_src, sizeof(uint8_t) * MAX_WEIGHTS_PER_BLOCK);
		}
	}

	
	// then iterate over the 4 believed-to-be-best modes to find out which one is
	// actually best.
	for (int l = 0; l < ewp.max_refinement_iters; l++)
	{
		for (int blk_idx = 0; blk_idx < batch_size; blk_idx++)
		{
			if (blk_stat[blk_idx] & skip_mode)
				continue;

			const imageblock * blk = &blk_batch[blk_idx];
			const error_weight_block *ewb = &ewb_batch[blk_idx];
			symbolic_compressed_block *scb = &scb_candidates[SCB_CANDIDATES * blk_idx];

			int partition_index = 0;
			if (partition_count > 1)
			{
				partition_index = partition_indices_2planes_batch[blk_idx * PARTITION_CANDIDATES + partition_offset];
				separate_component = partition_index >> PARTITION_BITS;
				partition_index = partition_index & (PARTITION_COUNT - 1);
			}
			const partition_info *pi = &ptab[partition_index];

			for (int i = 0; i < 4; i++)
			{
				if (tmpplanes.scb_stat[blk_idx] & (0x11 << i))
					continue;

				int scb_idx = SCB_CANDIDATES * blk_idx + i;

				const int weight_mode = tmpplanes.per_scb.weight_mode[scb_idx];
				const int decimation_mode = tmpplanes.per_scb.decimation_mode[scb_idx];
				const int color_quantization_level = tmpplanes.per_scb.color_quantization_level[scb_idx];
				const int color_quantization_level_mod = tmpplanes.per_scb.color_quantization_level_mod[scb_idx];
				const int *partition_format_specifiers = &tmpplanes.per_scb.partition_format_specifiers[MAX_PARTITIONS * scb_idx];
				const int weight_quantization_mode = sorted_bsd->block_modes[weight_mode].quantization_mode;
				const decimation_table *it = ixtab3[decimation_mode];

				const uint8_t *u8_weight1_src = &tmpplanes.per_scb.u8_qdq_weights[MAX_WEIGHTS_PER_BLOCK * scb_idx];
				const uint8_t *u8_weight2_src = &tmpplanes.per_scb.u8_qdq_weights[MAX_WEIGHTS_PER_BLOCK * scb_idx + MAX_WEIGHTS_PER_BLOCK / 2];

				// recompute the ideal color endpoints before storing them.
				endpoints *ep = &tmpplanes.per_scb.ep[scb_idx];
				float4 rgbs_colors[4];
				float4 rgbo_colors[4];
				float2 lum_intervals[4];


				recompute_ideal_colors(xdim, ydim, zdim, weight_quantization_mode, ep, rgbs_colors, rgbo_colors, lum_intervals, u8_weight1_src, u8_weight2_src, separate_component, pi, it, blk, ewb);

				// store the colors for the block
				for (int j = 0; j < partition_count; j++)
				{
					scb[i].color_formats[j] = pack_color_endpoints(decode_mode,
						ep->endpt0[j],
						ep->endpt1[j],
						rgbs_colors[j], rgbo_colors[j], lum_intervals[j], partition_format_specifiers[j], scb[i].color_values[j], color_quantization_level);
				}

				// if all the color endpoint modes are the same, we get a few more
				// bits to store colors; let's see if we can take advantage of this:
				// requantize all the colors and see if the endpoint modes remain the same;
				// if they do, then exploit it.
				scb[i].color_formats_matched = 0;

				if ((partition_count >= 2 && scb[i].color_formats[0] == scb[i].color_formats[1]
					&& color_quantization_level != color_quantization_level_mod)
					&& (partition_count == 2 || (scb[i].color_formats[0] == scb[i].color_formats[2] && (partition_count == 3 || (scb[i].color_formats[0] == scb[i].color_formats[3])))))
				{
					int colorvals[4][12];
					int color_formats_mod[4];
					for (int j = 0; j < partition_count; j++)
					{
						color_formats_mod[j] = pack_color_endpoints(decode_mode,
							ep->endpt0[j],
							ep->endpt1[j],
							rgbs_colors[j], rgbo_colors[j], lum_intervals[j], partition_format_specifiers[j], colorvals[j], color_quantization_level_mod);
					}
					if (color_formats_mod[0] == color_formats_mod[1]
						&& (partition_count == 2 || (color_formats_mod[0] == color_formats_mod[2] && (partition_count == 3 || (color_formats_mod[0] == color_formats_mod[3])))))
					{
						scb[i].color_formats_matched = 1;
						for (int j = 0; j < 4; j++)
							for (int k = 0; k < 12; k++)
								scb[i].color_values[j][k] = colorvals[j][k];
						for (int j = 0; j < 4; j++)
							scb[i].color_formats[j] = color_formats_mod[j];
					}
				}


				// store header fields
				scb[i].partition_count = partition_count;
				scb[i].partition_index = partition_index;
				scb[i].color_quantization_level = scb[i].color_formats_matched ? color_quantization_level_mod : color_quantization_level;
				scb[i].block_mode = sorted_bsd->block_modes[weight_mode].block_mode;
				scb[i].plane2_color_component = separate_component;
				scb[i].error_block = 0;

				if (scb[i].color_quantization_level < 4)
				{
					scb[i].error_block = 1;	// should never happen, but cannot prove it impossible
				}

			}
		}


		for (int scb_idx = 0; scb_idx < batch_size * SCB_CANDIDATES; scb_idx++)
		{
			int blk_idx = scb_idx >> 2;
			int scb_off = scb_idx & 0x03;
			if (tmpplanes.scb_stat[blk_idx] & (0x11 << scb_off))
				continue;
		
			const imageblock * blk = &blk_batch[blk_idx];
			const error_weight_block *ewb = &ewb_batch[blk_idx];

			uint8_t *u8_weight1_src = &tmpplanes.per_scb.u8_qdq_weights[MAX_WEIGHTS_PER_BLOCK * scb_idx];
			uint8_t *u8_weight2_src = &tmpplanes.per_scb.u8_qdq_weights[MAX_WEIGHTS_PER_BLOCK * scb_idx + MAX_WEIGHTS_PER_BLOCK / 2];

			int adjustments = realign_weights(decode_mode,
				xdim, ydim, zdim,
				blk, ewb, &scb_candidates[scb_idx],
				u8_weight1_src,
				u8_weight2_src);

			if (adjustments == 0)
				tmpplanes.scb_stat[blk_idx] |= 0x10 << scb_off;
		}
	}

	for (int scb_idx = 0; scb_idx < batch_size * SCB_CANDIDATES; scb_idx++)
	{
		int blk_idx = scb_idx >> 2;
		int scb_off = scb_idx & 0x03;
		if (tmpplanes.scb_stat[blk_idx] & (1 << scb_off))
			continue;

		int decimation_mode = tmpplanes.per_scb.decimation_mode[scb_idx];
		const decimation_table *it = ixtab3[decimation_mode];

		const uint8_t *u8_weight1_src = &tmpplanes.per_scb.u8_qdq_weights[MAX_WEIGHTS_PER_BLOCK * scb_idx];
		const uint8_t *u8_weight2_src = &tmpplanes.per_scb.u8_qdq_weights[MAX_WEIGHTS_PER_BLOCK * scb_idx + MAX_WEIGHTS_PER_BLOCK / 2];
		int weights_to_copy = it->num_weights;

		for (int j = 0; j < weights_to_copy; j++)
		{
			scb_candidates[scb_idx].plane1_weights[j] = u8_weight1_src[j];
			scb_candidates[scb_idx].plane2_weights[j] = u8_weight2_src[j];
		}
	}
}



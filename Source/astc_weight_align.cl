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
 *	@brief	Angular-sum algorithm for weight alignment.
 *
 *			This algorithm works as follows:
 *			* we compute a complex number P as (cos s*i, sin s*i) for each
 *			  weight, where i is the input value and s is a scaling factor
 *			  based on the spacing between the weights.
 *			* we then add together complex numbers for all the weights.
 *			* we then compute the length and angle of the resulting sum.
 *
 *			This should produce the following results:
 *			* perfect alignment results in a vector whose length is equal to
 *			  the sum of lengths of all inputs
 *			* even distribution results in a vector of length 0.
 *			* all samples identical results in perfect alignment for every
 *			  scaling.
 *
 *			For each scaling factor within a given set, we compute an alignment
 *			factor from 0 to 1. This should then result in some scalings standing
 *			out as having particularly good alignment factors; we can use this to
 *			produce a set of candidate scale/shift values for various quantization
 *			levels; we should then actually try them and see what happens.
 *
 *			Assuming N quantization steps, the scaling factor becomes s=2*PI*(N-1);
 *			we should probably have about 1 scaling factor for every 1/4
 *			quantization step (perhaps 1/8 for low levels of quantization)
 */ 
/*----------------------------------------------------------------------------*/ 

#include "astc_codec_internals_ocl.h"

constant float angular_steppings[] = {
	1.0, 1.125,
	1.25, 1.375,
	1.5, 1.625,
	1.75, 1.875,

	2.0, 2.25, 2.5, 2.75,
	3.0, 3.25, 3.5, 3.75,
	4.0, 4.25, 4.5, 4.75,
	5.0, 5.25, 5.5, 5.75,
	6.0, 6.25, 6.5, 6.75,
	7.0, 7.25, 7.5, 7.75,

	8.0, 8.5,
	9.0, 9.5,
	10.0, 10.5,
	11.0, 11.5,
	12.0, 12.5,
	13.0, 13.5,
	14.0, 14.5,
	15.0, 15.5,
	16.0, 16.5,
	17.0, 17.5,
	18.0, 18.5,
	19.0, 19.5,
	20.0, 20.5,
	21.0, 21.5,
	22.0, 22.5,
	23.0, 23.5,
	24.0, 24.5,
	25.0, 25.5,
	26.0, 26.5,
	27.0, 27.5,
	28.0, 28.5,
	29.0, 29.5,
	30.0, 30.5,
	31.0, 31.5,
	32.0, 32.5,
	33.0, 33.5,
	34.0, 34.5,
	35.0, 35.5,
};

#define ANGULAR_STEPS ((int)(sizeof(angular_steppings)/sizeof(angular_steppings[0])))

constant float stepsizes[ANGULAR_STEPS] = {
	1.000000000f, 0.888888896f, 0.800000012f, 0.727272749f, 0.666666687f, 0.615384638f, 0.571428597f, 0.533333361f,
	0.500000000f, 0.444444448f, 0.400000006f, 0.363636374f, 0.333333343f, 0.307692319f, 0.285714298f, 0.266666681f,
	0.250000000f, 0.235294119f, 0.222222224f, 0.210526317f, 0.200000003f, 0.190476194f, 0.181818187f, 0.173913047f,
	0.166666672f, 0.159999996f, 0.153846160f, 0.148148149f, 0.142857149f, 0.137931034f, 0.133333340f, 0.129032254f,
	0.125000000f, 0.117647059f, 0.111111112f, 0.105263159f, 0.100000001f, 0.0952380970f, 0.0909090936f, 0.0869565234f,
	0.0833333358f, 0.0799999982f, 0.0769230798f, 0.0740740746f, 0.0714285746f, 0.0689655170f, 0.0666666701f, 0.0645161271f,
	0.0625000000f, 0.0606060624f, 0.0588235296f, 0.0571428575f, 0.0555555560f, 0.0540540554f, 0.0526315793f, 0.0512820520f,
	0.0500000007f, 0.0487804860f, 0.0476190485f, 0.0465116277f, 0.0454545468f, 0.0444444455f, 0.0434782617f, 0.0425531901f,
	0.0416666679f, 0.0408163257f, 0.0399999991f, 0.0392156877f, 0.0384615399f, 0.0377358496f, 0.0370370373f, 0.0363636352f,
	0.0357142873f, 0.0350877196f, 0.0344827585f, 0.0338983051f, 0.0333333351f, 0.0327868834f, 0.0322580636f, 0.0317460336f,
	0.0312500000f, 0.0307692308f, 0.0303030312f, 0.0298507456f, 0.0294117648f, 0.0289855078f, 0.0285714287f, 0.0281690136f
};

constant int max_angular_steps_needed_for_quant_level[13] = { 8, 12, 16, 20, 24, 32, 36, 40, 48, 56, 64, 82, 87 };

#define SINCOS_STEPS 64


typedef union
{
	float f;
	int32_t s;
	uint32_t u;
} if32;


// function to compute angular sums; then, from the
// angular sums, compute alignment factor and offset.

/* static inline */
void compute_angular_offsets(int samplecount, global const float *samples, global const float *sample_weights, int max_angular_steps, float *offsets)
{
	int i, j;

	float anglesum_x[ANGULAR_STEPS];
	float anglesum_y[ANGULAR_STEPS];

	for (i = 0; i < max_angular_steps; i++)
	{
		anglesum_x[i] = 0;
		anglesum_y[i] = 0;
	}


	// compute the angle-sums.
	for (i = 0; i < samplecount; i++)
	{
#ifdef MIMIC_HOST
		float sample = samples[i];
		float sample_weight = sample_weights[i];
		if32 p;
		p.f = (sample * (SINCOS_STEPS - 1.0f)) + 12582912.0f;
		unsigned int isample = p.u & 0x3F;

		for (j = 0; j < max_angular_steps; j++)
		{
			float cp = cos((2.0f * M_PI / (SINCOS_STEPS - 1.0f)) * angular_steppings[j] * isample);
			float sp = sin((2.0f * M_PI / (SINCOS_STEPS - 1.0f)) * angular_steppings[j] * isample);

			anglesum_x[j] += cp * sample_weight;
			anglesum_y[j] += sp * sample_weight;
		}
#else
		float sample = samples[i];
		float sample_weight = sample_weights[i];
		
		for (j = 0; j < max_angular_steps; j++)
		{
			float angle = (2.0f * M_PI ) * angular_steppings[j] * sample;
			float cp = native_cos(angle);
			float sp = native_sin(angle);

			anglesum_x[j] += cp * sample_weight;
			anglesum_y[j] += sp * sample_weight;
		}
#endif // MIMIC_HOST
	}

	// postprocess the angle-sums
	for (i = 0; i < max_angular_steps; i++)
	{
		float angle = atan2(anglesum_y[i], anglesum_x[i]);	// positive angle -> positive offset
		offsets[i] = angle * (stepsizes[i] * (1.0f / (2.0f * (float)M_PI)));
	}
}



// for a given step-size and a given offset, compute the
// lowest and highest weight that results from quantizing using the stepsize & offset.
// also, compute the resulting error.


/* static inline */
void compute_lowest_and_highest_weight(int samplecount, global const float *samples, global const float *sample_weights,
									  int max_angular_steps, const float *offsets,
									  int8_t * lowest_weight, int8_t * highest_weight,
									  float *error, float *cut_low_weight_error, float *cut_high_weight_error)
{
	int i;
	int sp;

	for (sp = 0; sp < max_angular_steps; sp++)
	{
		float rcp_stepsize = angular_steppings[sp];
		float offset = offsets[sp];

		float scaled_offset = rcp_stepsize * offset;


		float wt = sample_weights[0];
		if32 p;
		float sval = (samples[0] * rcp_stepsize) - scaled_offset;
		p.f = sval + 12582912.0f;	// FP representation abuse to avoid floor() and float->int conversion. TODO: replace with float rounding
		float isval = p.f - 12582912.0f;
		float dif = sval - isval;

		float errval = (dif * wt) * dif;

		int sh = rint(sval);
		int idx_bias12 = clamp(sh, -12, 43); // TODO: is clamp needed?

		int minidx_bias12 = idx_bias12;
		int maxidx_bias12 = idx_bias12;

		float error_cut_low = (1.0f - 2.0f * dif) * wt;
		float error_cut_high = (1.0f + 2.0f * dif) * wt;

		for (i = 1; i < samplecount; i++)
		{
			float wt = sample_weights[i];
			if32 p;
			float sval = (samples[i] * rcp_stepsize) - scaled_offset;
			p.f = sval + 12582912.0f;	// FP representation abuse to avoid floor() and float->int conversion. TODO: replace with float rounding
			float isval = p.f - 12582912.0f;
			float dif = sval - isval;

			errval += (dif * wt) * dif;

			int sh = rint(sval);
			int idx_bias12 = clamp(sh, -12, 43); // TODO: is clamp needed?

			//zeroing min/max errors if min/max index have changed
			error_cut_low *= (idx_bias12 >= minidx_bias12);
			error_cut_high *= (idx_bias12 <= maxidx_bias12);

			minidx_bias12 = min(minidx_bias12, idx_bias12);
			maxidx_bias12 = max(maxidx_bias12, idx_bias12);

			error_cut_low += (1.0f - 2.0f * dif) * wt * (idx_bias12 == minidx_bias12);
			error_cut_high += (1.0f + 2.0f * dif) * wt * (idx_bias12 == maxidx_bias12);
		}


		lowest_weight[sp] = minidx_bias12;
		highest_weight[sp] = maxidx_bias12;
		error[sp] = errval;

		// the cut_(lowest/highest)_weight_error indicate the error that results from
		// forcing samples that should have had the (lowest/highest) weight value
		// one step (up/down).
		cut_low_weight_error[sp] = error_cut_low;
		cut_high_weight_error[sp] = error_cut_high;
	}


	for (sp = 0; sp < max_angular_steps; sp++)
	{
		float errscale = stepsizes[sp] * stepsizes[sp];
		error[sp] *= errscale;
		cut_low_weight_error[sp] *= errscale;
		cut_high_weight_error[sp] *= errscale;
	}
}



// main function for running the angular algorithm.


void compute_angular_endpoints_for_quantization_levels(int samplecount, global const float *samples, global const float *sample_weights, int max_quantization_level, float low_value[12], float high_value[12])
{
	int i;


	max_quantization_level++;	// Temporarily increase level - needs refinement

	const int quantization_steps_for_level[13] = { 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 33, 36 };
	int max_quantization_steps = quantization_steps_for_level[max_quantization_level];

	float offsets[ANGULAR_STEPS];

	int max_angular_steps = max_angular_steps_needed_for_quant_level[max_quantization_level];

	compute_angular_offsets(samplecount, samples, sample_weights, max_angular_steps, offsets);


	int8_t lowest_weight[ANGULAR_STEPS];
	int8_t highest_weight[ANGULAR_STEPS];
	float error[ANGULAR_STEPS];

	float cut_low_weight_error[ANGULAR_STEPS];
	float cut_high_weight_error[ANGULAR_STEPS];

	compute_lowest_and_highest_weight(samplecount, samples, sample_weights, max_angular_steps, offsets, lowest_weight, highest_weight, error, cut_low_weight_error, cut_high_weight_error);


	// for each quantization level, find the best error terms.
	float best_errors[40];
	int best_scale[40];
	uint8_t cut_low_weight[40];

	for (i = 0; i < (max_quantization_steps + 4); i++)
	{
		best_errors[i] = 1e30f;
		best_scale[i] = -1;	// Indicates no solution found
		cut_low_weight[i] = 0;
	}



	for (i = 0; i < max_angular_steps; i++)
	{
		int samplecount = highest_weight[i] - lowest_weight[i] + 1;
		if (samplecount >= (max_quantization_steps + 4))
		{
			continue;
		}
		if (samplecount < 2)
			samplecount = 2;

		if (best_errors[samplecount] > error[i])
		{
			best_errors[samplecount] = error[i];
			best_scale[samplecount] = i;
			cut_low_weight[samplecount] = 0;
		}

		float error_cut_low = error[i] + cut_low_weight_error[i];
		float error_cut_high = error[i] + cut_high_weight_error[i];
		float error_cut_low_high = error[i] + cut_low_weight_error[i] + cut_high_weight_error[i];

		if (best_errors[samplecount - 1] > error_cut_low)
		{
			best_errors[samplecount - 1] = error_cut_low;
			best_scale[samplecount - 1] = i;
			cut_low_weight[samplecount - 1] = 1;
		}

		if (best_errors[samplecount - 1] > error_cut_high)
		{
			best_errors[samplecount - 1] = error_cut_high;
			best_scale[samplecount - 1] = i;
			cut_low_weight[samplecount - 1] = 0;
		}

		if (best_errors[samplecount - 2] > error_cut_low_high)
		{
			best_errors[samplecount - 2] = error_cut_low_high;
			best_scale[samplecount - 2] = i;
			cut_low_weight[samplecount - 2] = 1;
		}

	}

	// if we got a better error-value for a low samplecount than for a high one,
	// use the low-samplecount error value for the higher samplecount as well.
	for (i = 3; i <= max_quantization_steps; i++)
	{
		if (best_errors[i] > best_errors[i - 1])
		{
			best_errors[i] = best_errors[i - 1];
			best_scale[i] = best_scale[i - 1];
			cut_low_weight[i] = cut_low_weight[i - 1];
		}
	}


	max_quantization_level--;	// Decrease level again (see corresponding ++, above)

	//const int ql_weights[12] = { 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 33 };
	for (i = 0; i <= max_quantization_level; i++)
	{
		int q = quantization_steps_for_level[i];//ql_weights[i];
		int bsi = best_scale[q];

		// Did we find anything?
		if(bsi < 0)
		{
			//printf("ERROR: Unable to find an encoding within the specified error limits. Please revise the error limit values and try again.\n");
			//exit(1);
		}

		float stepsize = stepsizes[bsi];
		int lwi = lowest_weight[bsi] + cut_low_weight[q];
		int hwi = lwi + q - 1;
		float offset = offsets[bsi];

		low_value[i] = offset + lwi * stepsize;
		high_value[i] = offset + hwi * stepsize;
	}

}


// helper functions that will compute ideal angular-endpoints
// for a given set of weights and a given block size descriptors

__kernel
void compute_angular_endpoints_1plane(global const uint8_t *blk_stat, global const block_size_descriptor_sorted * bsds,
									  global const float *g_decimated_quantized_weights, global const float *g_decimated_weights,
									  global float g_low_value[WLIMIT_1PLANE], global float g_high_value[WLIMIT_1PLANE])
{
	uint blk_idx = get_global_id(0);
	if (blk_stat[blk_idx] & BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF)
		return;
	
	global const float* decimated_quantized_weights = g_decimated_quantized_weights + DLIMIT_1PLANE * MAX_WEIGHTS_PER_BLOCK * blk_idx;
	global const float* decimated_weights = g_decimated_weights + DLIMIT_1PLANE * MAX_WEIGHTS_PER_BLOCK * blk_idx;
	global float* low_value = g_low_value + WLIMIT_1PLANE * blk_idx;
	global float* high_value = g_high_value + WLIMIT_1PLANE * blk_idx;

	int i;
	float low_values[DLIMIT_1PLANE][12];
	float high_values[DLIMIT_1PLANE][12];

	for (i = 0; i < DLIMIT_1PLANE; i++) // DLIMIT_1PLANE = ewp->decimation_mode_limit_1plane
	{
		int samplecount = bsds->decimation_mode_samples[i];
		int quant_mode = bsds->decimation_mode_maxprec[i];

		compute_angular_endpoints_for_quantization_levels(samplecount,
														  decimated_quantized_weights + i * MAX_WEIGHTS_PER_BLOCK,
														  decimated_weights + i * MAX_WEIGHTS_PER_BLOCK, quant_mode, low_values[i], high_values[i]);
	}

	for (i = 0; i < WLIMIT_1PLANE; i++) // WLIMIT_1PLANE = ewp->weight_mode_limit_1plane
	{
		int quant_mode = bsds->block_modes[i].quantization_mode;
		int decim_mode = bsds->block_modes[i].sorted_decimation_mode;

		low_value[i] = low_values[decim_mode][quant_mode];
		high_value[i] = high_values[decim_mode][quant_mode];
	}

}


__kernel
void compute_angular_endpoints_2planes(global const uint8_t *blk_stat, global const block_size_descriptor_sorted * bsds,
									   global const float *g_decimated_quantized_weights,
									   global const float *g_decimated_weights,
									   global float g_low_value1[WLIMIT_2PLANES], global float g_high_value1[WLIMIT_2PLANES], global float g_low_value2[WLIMIT_2PLANES], global float g_high_value2[WLIMIT_2PLANES])
{
	uint blk_idx = get_global_id(0);
	if (blk_stat[blk_idx] & BLOCK_STAT_TEXEL_AVG_ERROR_CUTOFF)
		return;

	global const float* decimated_quantized_weights = g_decimated_quantized_weights + DLIMIT_2PLANES * MAX_WEIGHTS_PER_BLOCK * blk_idx;
	global const float* decimated_weights = g_decimated_weights + DLIMIT_2PLANES * MAX_WEIGHTS_PER_BLOCK * blk_idx;
	global float* low_value1 = g_low_value1 + WLIMIT_2PLANES * blk_idx;
	global float* low_value2 = g_low_value2 + WLIMIT_2PLANES * blk_idx;
	global float* high_value1 = g_high_value1 + WLIMIT_2PLANES * blk_idx;
	global float* high_value2 = g_high_value2 + WLIMIT_2PLANES * blk_idx;

	int i;
	float low_values1[DLIMIT_2PLANES][12];
	float high_values1[DLIMIT_2PLANES][12];
	float low_values2[DLIMIT_2PLANES][12];
	float high_values2[DLIMIT_2PLANES][12];

	for (i = 0; i < DLIMIT_2PLANES; i++) // DLIMIT_2PLANES = ewp->decimation_mode_limit_2planes
	{
		int samplecount = bsds->decimation_mode_samples[i];
		int quant_mode = bsds->decimation_mode_maxprec[i];

		compute_angular_endpoints_for_quantization_levels(samplecount,
														  decimated_quantized_weights + 2 * i * (MAX_WEIGHTS_PER_BLOCK / 2),
														  decimated_weights + 2 * i * (MAX_WEIGHTS_PER_BLOCK / 2), quant_mode, low_values1[i], high_values1[i]);

		compute_angular_endpoints_for_quantization_levels(samplecount,
														  decimated_quantized_weights + (2 * i + 1) * (MAX_WEIGHTS_PER_BLOCK / 2),
														  decimated_weights + (2 * i + 1) * (MAX_WEIGHTS_PER_BLOCK / 2), quant_mode, low_values2[i], high_values2[i]);

	}

	for (i = 0; i < WLIMIT_2PLANES; i++) // WLIMIT_2PLANES = ewp->weight_mode_limit_2planes
	{
		int quant_mode = bsds->block_modes[i].quantization_mode;
		int decim_mode = bsds->block_modes[i].sorted_decimation_mode;

		low_value1[i] = low_values1[decim_mode][quant_mode];
		high_value1[i] = high_values1[decim_mode][quant_mode];
		low_value2[i] = low_values2[decim_mode][quant_mode];
		high_value2[i] = high_values2[decim_mode][quant_mode];
	}
}

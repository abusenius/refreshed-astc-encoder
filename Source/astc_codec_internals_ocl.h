/*----------------------------------------------------------------------------*/  
/**
 *
 *	@brief	Internal function and data declarations for ASTC codec.
 *          Include file for OpenCL C kernels
 */ 
/*----------------------------------------------------------------------------*/ 

#ifndef ASTC_CODEC_INTERNALS_OCL_INCLUDED
#define ASTC_CODEC_INTERNALS_OCL_INCLUDED

#ifdef CL_VERSION_1_2
    // included from OpenCL kernel
	#define OPENCL_C_KERNEL

	typedef signed char int8_t;
	typedef signed short int16_t;
	typedef ulong uint64_t;
	typedef unsigned char uint8_t;
	typedef unsigned short uint16_t;
#else
	// just to make IDE happy with OpenCL keywords and macros
	static_assert(false, "This file could be included only from OpenCL kernel");
	#define __global
	#define __local
	#define __kernel
    #define get_global_id(i) i
	#define sqrt(i) i

	#define XDIM 4
	#define YDIM 4
	#define ZDIM 1
	#define TEXELS_PER_BLOCK (XDIM * YDIM * ZDIM)
	#define WEIGHT_IMPRECISION_ESTIM_SQUARED 0.03f
#endif

#include "astc_codec_internals.h"

//declarations from mathlib.h

// parametric line, 2D: The line is given by line = a + b*t.
typedef struct
{
	float2 a;
	float2 b;
} line2;

// paramtric line, 3D
typedef struct
{
	float3 a;
	float3 b;
} line3;

typedef struct
{
	float4 a;
	float4 b;
} line4;

// plane/hyperplane defined by a point and a normal vector
typedef struct
{
	float3 root_point;
	float3 normal;				// normalized
} plane_3d;

typedef struct
{
	float4 root_point;
	float4 normal;				// normalized
} hyperplane_4d;



void compute_averages_and_directions_rgba(__global const partition_info * pt,
	__global const imageblock * blk,
	__global const error_weight_block * ewb,
	const float4 * color_scalefactors,
	float4 * averages, float4 * directions_rgba, float3 * directions_gba, float3 * directions_rba, float3 * directions_rga, float3 * directions_rgb);

void compute_averages_and_directions_rgb(__global const partition_info * pt,
	__global const imageblock * blk,
	__global const error_weight_block * ewb,
	const float4 * color_scalefactors, float3 * averages, float3 * directions_rgb, float2 * directions_rg, float2 * directions_rb, float2 * directions_gb);


float compute_error_squared_gb(__global const partition_info * pt, __global const imageblock * blk, __global const error_weight_block * ewb, const processed_line2 * plines, float *length_of_lines);
float compute_error_squared_rb(__global const partition_info * pt, __global const imageblock * blk, __global const error_weight_block * ewb, const processed_line2 * plines, float *length_of_lines);
float compute_error_squared_rg(__global const partition_info * pt, __global const imageblock * blk, __global const error_weight_block * ewb, const processed_line2 * plines, float *length_of_lines);
float compute_error_squared_gba(__global const partition_info * pt, __global const imageblock * blk, __global const error_weight_block * ewb, const processed_line3 * plines, float *length_of_lines);
float compute_error_squared_rba(__global const partition_info * pt, __global const imageblock * blk, __global const error_weight_block * ewb, const processed_line3 * plines, float *length_of_lines);
float compute_error_squared_rga(__global const partition_info * pt, __global const imageblock * blk, __global const error_weight_block * ewb, const processed_line3 * plines, float *length_of_lines);
float compute_error_squared_rgb(__global const partition_info * pt, __global const imageblock * blk, __global const error_weight_block * ewb, const processed_line3 * plines, float *length_of_lines);
float compute_error_squared_rgba(__global const partition_info * pt, __global const imageblock * blk, __global const error_weight_block * ewb, const processed_line4 * plines, float *length_of_lines);

#endif
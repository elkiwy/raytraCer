#ifndef __UTILS_H_
#define __UTILS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <OpenCL/cl.h>


typedef struct float3{
    float x;
    float y;
    float z;
}float3;




/**
 * Randoms
 * */

float random_double_unit();
float random_double(float min, float max);
cl_float3 random_in_unit_disk();
cl_float3 random_in_unit_sphere();



/**
 * Vector Math
 * */

float float3_length(cl_float3 v);
cl_float3 float3_normalize(cl_float3 v);
cl_float3 float3_cross(cl_float3 v1, cl_float3 v2);



/**
 * Perlin
 * */

#define PERLIN_POINT_COUNT 256
void perlin_init(cl_int3* perm_xyz, cl_float3* ranvec);



/**
 * OpenCL
 * */

cl_device_id create_device();
cl_program build_program(cl_context ctx, cl_device_id dev);





#endif // __UTILS_H_

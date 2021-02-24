#ifndef __UTILS_H_
#define __UTILS_H_

#include <math.h>
#include <stdlib.h>
#include <OpenCL/cl.h>

/**
 * Randoms
 * */

float random_double_unit() {
    return rand() / (RAND_MAX + 1.0);
}



float random_double(float min, float max) {
    return min + (max-min)*random_double_unit();
}



cl_float3 random_in_unit_disk() {
    while (1) {
        cl_float3 p = (cl_float3){{random_double(-1,1), random_double(-1,1), 0}};
        if((p.s[0]*p.s[0] + p.s[1]*p.s[1] + p.s[2]*p.s[2]) < 1){return p;}
    }
}



cl_float3 random_in_unit_sphere(){
    while (1) {
        cl_float3 p = {{random_double(-1,1), random_double(-1,1), random_double(-1,1)}};
        if((p.s[0]*p.s[0] + p.s[1]*p.s[1] + p.s[2]*p.s[2]) < 1){return p;}
    }
}





/**
 * Vector Math
 * */

float float3_length(cl_float3 v){
    return sqrt((v.s[0] * v.s[0]) + (v.s[1] * v.s[1]) + (v.s[2] * v.s[2]));
}



cl_float3 float3_normalize(cl_float3 v){
    float l = float3_length(v);
    return (cl_float3){{v.s[0]/l, v.s[1]/l, v.s[2]/l}};
}



cl_float3 float3_cross(cl_float3 v1, cl_float3 v2){
    return (cl_float3){{
        v1.s[1] * v2.s[2] - v1.s[2] * v2.s[1],
        v1.s[2] * v2.s[0] - v1.s[0] * v2.s[2],
        v1.s[0] * v2.s[1] - v1.s[1] * v2.s[0]
    }};
}








/**
 * Perlin
 * */

#define PERLIN_POINT_COUNT 256

void perlin_init(cl_int3* perm_xyz, cl_float3* ranvec){
    for (int i=0;i<PERLIN_POINT_COUNT; ++i){
        ranvec[i] = float3_normalize((cl_float3){{random_double(-1, 1), random_double(-1, 1), random_double(-1, 1)}});
    }
    for(int i=0;i<PERLIN_POINT_COUNT;++i){
        perm_xyz[i].s[0] = i;
        perm_xyz[i].s[1] = i;
        perm_xyz[i].s[2] = i;
    }
    for(int i=0;i<PERLIN_POINT_COUNT;++i){
        for(int j=0;j<3;j++){
            int target = rand() % PERLIN_POINT_COUNT;
            int tmp = perm_xyz[i].s[j];
            perm_xyz[i].s[j] = perm_xyz[target].s[j];
            perm_xyz[target].s[j] = tmp;
        }
    }
}








#endif // __UTILS_H_

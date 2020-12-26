#ifndef __VEC3_H_
#define __VEC3_H_

#include <stdio.h>

#include <math.h>

/**
 *
 * Vec3
 *
 * */

//Struct
typedef struct{double x; double y; double z;}vec3;
inline static void vec3_print(vec3* v){printf("Vec3: %f %f %f\n", v->x, v->y, v->z);}
inline static vec3 vec3_copy(vec3* v){vec3 r = {v->x, v->y, v->z};return r;}

//Vector operations without copies
inline static vec3 vec3_sum(vec3* v1, vec3* v2){vec3 result = {v1->x + v2->x, v1->y + v2->y, v1->z + v2->z}; return result;}
inline static vec3 vec3_sub(vec3* v1, vec3* v2){vec3 result = {v1->x - v2->x, v1->y - v2->y, v1->z - v2->z}; return result;}
inline static vec3 vec3_mul(vec3* v1, vec3* v2){vec3 result = {v1->x * v2->x, v1->y * v2->y, v1->z * v2->z}; return result;}
inline static vec3 vec3_mul_k(vec3* v, double k){vec3 result = {v->x * k, v->y * k, v->z * k}; return result;}
inline static vec3 vec3_div_k(vec3* v, double k){vec3 result = {v->x * 1/k, v->y * 1/k, v->z * 1/k}; return result;}
inline static double vec3_length(vec3* v){return sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));}
inline static double vec3_length_squared(vec3* v){return (v->x * v->x) + (v->y * v->y) + (v->z * v->z);}
inline static double vec3_dot(vec3* v1, vec3* v2){return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;}
inline static vec3 vec3_cross(vec3* v1, vec3* v2){vec3 result = {v1->y * v2->z - v1->z * v2->y, v1->z * v2->x - v1->x * v2->z, v1->x * v2->y - v1->y * v2->x}; return result;}
inline static vec3 vec3_unit(vec3* v1){vec3 res = vec3_div_k(v1, vec3_length(v1));return res;}
inline static vec3 vec3_flip(vec3* v1){vec3 res = {-v1->x, -v1->y, -v1->z};return res;}

//Vector operation with copies
inline static vec3 vec3c_sum(vec3 v1, vec3 v2){vec3 result = {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z}; return result;}
inline static vec3 vec3c_sub(vec3 v1, vec3 v2){vec3 result = {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z}; return result;}
inline static vec3 vec3c_mul(vec3 v1, vec3 v2){vec3 result = {v1.x * v2.x, v1.y * v2.y, v1.z * v2.z}; return result;}
inline static vec3 vec3c_mul_k(vec3 v, double k){vec3 result = {v.x * k, v.y * k, v.z * k}; return result;}
inline static vec3 vec3c_div_k(vec3 v, double k){vec3 result = {v.x * 1/k, v.y * 1/k, v.z * 1/k}; return result;}
inline static double vec3c_length(vec3 v){return sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));}
inline static double vec3c_length_squared(vec3 v){return (v.x * v.x) + (v.y * v.y) + (v.z * v.z);}
inline static double vec3c_dot(vec3 v1, vec3 v2){return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;}
inline static vec3 vec3c_cross(vec3 v1, vec3 v2){vec3 result = {v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x}; return result;}
inline static vec3 vec3c_unit(vec3 v1){vec3 res = vec3c_div_k(v1, vec3c_length(v1));return res;}
inline static vec3 vec3c_flip(vec3 v1){vec3 res = {-v1.x, -v1.y, -v1.z};return res;}

//Names
typedef vec3 point3;
typedef vec3 color;

#endif // __VEC3_H_

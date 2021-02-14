#ifndef __VEC3_H_
#define __VEC3_H_

#include <stdio.h>

#include <math.h>

#include "utils.h"

/**
 *
 * Vec3
 *
 * */

//Struct
typedef struct{double x; double y; double z;}vec3;
inline static void vec3_print(vec3* v){printf("Vec3: %.2f %.2f %.2f\n", v->x, v->y, v->z);}
inline static vec3 vec3_copy(vec3* v){vec3 r = {v->x, v->y, v->z};return r;}

//Vector operations without copies
inline static vec3 vec3_sum(vec3* v1, vec3* v2){return (vec3){v1->x + v2->x, v1->y + v2->y, v1->z + v2->z}; }
inline static vec3 vec3_sub(vec3* v1, vec3* v2){return (vec3){v1->x - v2->x, v1->y - v2->y, v1->z - v2->z}; }
inline static vec3 vec3_mul(vec3* v1, vec3* v2){return (vec3){v1->x * v2->x, v1->y * v2->y, v1->z * v2->z}; }
inline static vec3 vec3_mul_k(vec3* v, double k){return (vec3){v->x * k, v->y * k, v->z * k}; }
inline static vec3 vec3_div_k(vec3* v, double k){return (vec3){v->x * 1/k, v->y * 1/k, v->z * 1/k}; }
inline static double vec3_length(vec3* v){return sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));}
inline static double vec3_length_squared(vec3* v){return (v->x * v->x) + (v->y * v->y) + (v->z * v->z);}
inline static double vec3_dot(vec3* v1, vec3* v2){return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;}
inline static vec3 vec3_cross(vec3* v1, vec3* v2){return (vec3){v1->y * v2->z - v1->z * v2->y, v1->z * v2->x - v1->x * v2->z, v1->x * v2->y - v1->y * v2->x}; }
inline static vec3 vec3_unit(vec3* v1){return vec3_div_k(v1, vec3_length(v1));}
inline static vec3 vec3_flip(vec3* v1){return (vec3){-v1->x, -v1->y, -v1->z};}
inline static vec3 vec3_reflect(vec3* v, vec3* n){double a = vec3_dot(v, n) * 2.0; vec3 b = vec3_mul_k(n, a); return vec3_sub(v, &b);}

//Vector operation with copies
inline static vec3 vec3c_sum(vec3 v1, vec3 v2){return (vec3){v1.x + v2.x, v1.y + v2.y, v1.z + v2.z}; }
inline static vec3 vec3c_sub(vec3 v1, vec3 v2){return (vec3){v1.x - v2.x, v1.y - v2.y, v1.z - v2.z}; }
inline static vec3 vec3c_mul(vec3 v1, vec3 v2){return (vec3){v1.x * v2.x, v1.y * v2.y, v1.z * v2.z}; }
inline static vec3 vec3c_mul_k(vec3 v, double k){return (vec3){v.x * k, v.y * k, v.z * k}; }
inline static vec3 vec3c_div_k(vec3 v, double k){return (vec3){v.x * 1/k, v.y * 1/k, v.z * 1/k}; }
inline static double vec3c_length(vec3 v){return sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));}
inline static double vec3c_length_squared(vec3 v){return (v.x * v.x) + (v.y * v.y) + (v.z * v.z);}
inline static double vec3c_dot(vec3 v1, vec3 v2){return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;}
inline static vec3 vec3c_cross(vec3 v1, vec3 v2){return (vec3){v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x}; }
inline static vec3 vec3c_unit(vec3 v1){return vec3c_div_k(v1, vec3c_length(v1));}
inline static vec3 vec3c_flip(vec3 v1){return (vec3){-v1.x, -v1.y, -v1.z};}
inline static vec3 vec3c_reflect(vec3 v, vec3 n){double a = vec3c_dot(v, n) * 2.0; return vec3c_sub(v, vec3c_mul_k(n, a));}


///Refract a vector with a normal and a value etai_over_etat
inline static vec3 vec3c_refract(vec3 uv, vec3 n, double etai_over_etat){
    double cos_theta = fmin(vec3c_dot(vec3_flip(&uv), n), 1.0);
    vec3 r_out_perp =  vec3c_mul_k(vec3c_sum(uv, vec3c_mul_k(n, cos_theta)), etai_over_etat);
    double k =  -sqrt(fabs(1.0 - vec3_length_squared(&r_out_perp)));
    vec3 r_out_parallel =  vec3_mul_k(&n, k);
    return vec3c_sum(r_out_perp, r_out_parallel);
}

///Check if a vector is near zero
int vec3_is_near_zero(vec3* v){return (fabs(v->x) < 1e-8) && (fabs(v->y) < 1e-8) && (fabs(v->z) < 1e-8);}

///Create a random vector of 0-1 components
inline static vec3 vec3_random(){
    return (vec3){random_double(), random_double(), random_double()};
}

///Create a random vector of min-max components
inline static vec3 vec3_random_scaled(double min, double max){
    return (vec3){random_double_scaled(min, max), random_double_scaled(min, max), random_double_scaled(min, max)};
}

///Create a random vector inside a unit sphere
inline static vec3 vec3_random_in_unit_sphere(){
    while(1){
        vec3 p = vec3_random_scaled(-1, 1);
        if (vec3_length_squared(&p) >= 1) continue;
        return p;
    }
}

///Random a random UNIT vector inside a unit sphere
inline static vec3 vec3_random_unit_vector_in_unit_sphere(){
    return vec3c_unit(vec3_random_in_unit_sphere());
}

///Create a random unit vector in hemisphere
vec3 vec3_random_in_hemisphere(vec3* normal){
    vec3 in_unit_sphere = vec3_random_in_unit_sphere();
    if (vec3_dot(&in_unit_sphere, normal) > 0.0){return in_unit_sphere;
    }else{return vec3_flip(&in_unit_sphere);}
}

///Random vector in unit disk
vec3 random_in_unit_disk(){
    while(1){
        vec3 p = {random_double_scaled(-1, 1), random_double_scaled(-1, 1), 0};
        if(vec3_length_squared(&p) >= 1) continue;
        return p;
    }
}





//Names
typedef vec3 point3;
typedef vec3 color;

#endif // __VEC3_H_

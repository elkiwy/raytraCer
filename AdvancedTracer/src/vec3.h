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
void vec3_print(vec3* v);
vec3 vec3_copy(vec3* v);

//Vector operations without copies
vec3 vec3_sum(vec3* v1, vec3* v2);
vec3 vec3_sub(vec3* v1, vec3* v2);
vec3 vec3_mul(vec3* v1, vec3* v2);
vec3 vec3_mul_k(vec3* v, double k);
vec3 vec3_div_k(vec3* v, double k);
double vec3_length(vec3* v);
double vec3_length_squared(vec3* v);
double vec3_dot(vec3* v1, vec3* v2);
vec3 vec3_cross(vec3* v1, vec3* v2);
vec3 vec3_unit(vec3* v1);
vec3 vec3_flip(vec3* v1);
vec3 vec3_reflect(vec3* v, vec3* n);


//Vector operation with copies
vec3 vec3c_sum(vec3 v1, vec3 v2);
vec3 vec3c_sub(vec3 v1, vec3 v2);
vec3 vec3c_mul(vec3 v1, vec3 v2);
vec3 vec3c_mul_k(vec3 v, double k);
vec3 vec3c_div_k(vec3 v, double k);
double vec3c_length(vec3 v);
double vec3c_length_squared(vec3 v);
double vec3c_dot(vec3 v1, vec3 v2);
vec3 vec3c_cross(vec3 v1, vec3 v2);
vec3 vec3c_unit(vec3 v1);
vec3 vec3c_flip(vec3 v1);
vec3 vec3c_reflect(vec3 v, vec3 n);


///Refract 
vec3 vec3c_refract(vec3 uv, vec3 n, double etai_over_etat);



///Check if a vector is near zero
int vec3_is_near_zero(vec3* v);

///Create a random vector of 0-1 components
vec3 vec3_random();

///Create a random vector of min-max components
vec3 vec3_random_scaled(double min, double max);

///Create a random vector inside a unit sphere
vec3 vec3_random_in_unit_sphere();

///Random a random UNIT vector inside a unit sphere
vec3 vec3_random_unit_vector_in_unit_sphere();

///Create a random unit vector in hemisphere
vec3 vec3_random_in_hemisphere(vec3* normal);

///Random vector in unit disk
vec3 random_in_unit_disk();





//Names
typedef vec3 point3;
typedef vec3 color;

#endif // __VEC3_H_

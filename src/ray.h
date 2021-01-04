#ifndef __RAY_H_
#define __RAY_H_


#include <stdio.h>

#include "vec3.h"


/**
 *
 * Ray
 *
 * */

///Ray structure
typedef struct{point3 orig; vec3 dir; double time;} ray;

///Move the origin point through the ray for T distance
point3 ray_at(ray* r, double t);

///Print a ray object
void ray_print(char* s, ray* r);


#endif // __RAY_H_

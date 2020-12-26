#ifndef __RAY_H_
#define __RAY_H_


#include <stdio.h>

#include "vec3.h"


/**
 *
 * Ray
 *
 * */

typedef struct{point3 orig; vec3 dir;} ray;

inline static point3 ray_at(ray* r, double t){
    vec3 a = vec3_mul_k(&r->dir, t);
    return vec3_sum(&r->orig, &a);
}

void ray_print(ray* r){
    printf("Ray: orig: [%f %f %f] dir: [%f %f %f]", r->orig.x, r->orig.y, r->orig.z, r->dir.x, r->dir.y, r->dir.z);
}


#endif // __RAY_H_

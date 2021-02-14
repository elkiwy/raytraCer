#ifndef __PERLIN_H_
#define __PERLIN_H_

#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "vec3.h"

#define POINT_COUNT 256


typedef struct{
    //There are all arrays of POINT_COUNT, but we only store the point in the struct for efficency
    vec3* ranvec; int* perm_x; int* perm_y; int* perm_z;
}perlin;



perlin* perlin_init();
void perlin_free(perlin* p);
double perlin_noise(perlin* p, point3 point);
double turb_rec(perlin* perl, point3 p, int depth);
double perlin_turb(perlin* perl, point3 p);

#endif // __PERLIN_H_

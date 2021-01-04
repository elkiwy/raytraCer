#ifndef __AABB_H_
#define __AABB_H_


#include "utils.h"
#include "vec3.h"
#include "ray.h"



typedef struct{
    point3 minimum;
    point3 maximum;
}aabb;


aabb* aabb_init(point3 a, point3 b);
void aabb_free(aabb* b);
void aabb_print(aabb* b);
int aabb_hit(aabb* box, ray* r, double tmin, double tmax);
aabb surrounding_box(aabb box0, aabb box1);




#endif // __AABB_H_

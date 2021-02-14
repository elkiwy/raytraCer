#ifndef __AABB_H_
#define __AABB_H_


#include "utils.h"
#include "vec3.h"
#include "ray.h"

/**
* Axis Aligned Bounding Box class
*  used for BVH Tree optimization and space object grouping
*/

///AABB struct
typedef struct{point3 minimum; point3 maximum;}aabb;


//Inits
aabb* aabb_init(point3 a, point3 b);
aabb surrounding_box(aabb box0, aabb box1);

//Deinit
void aabb_free(aabb* b);

//Prints
void aabb_print(aabb* b);

//Features
int aabb_hit(aabb* box, ray* r, double tmin, double tmax);




#endif // __AABB_H_

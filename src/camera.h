#ifndef __CAMERA_H_
#define __CAMERA_H_

#include <stdlib.h>
#include "vec3.h"
#include "ray.h"

/**
 *
 * Camera
 *
 * */

typedef struct {
    point3 origin;
    point3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
}camera;

camera* camera_new(){
    camera* c = malloc(sizeof(camera));

    const double ASPECT_RATIO = 16.0/9.0;
    double viewport_height = 2.0;
    double viewport_width = ASPECT_RATIO * viewport_height;
    double focal_length = 1.0;

    c->origin = (vec3){0, 0, 0};
    c->horizontal = (vec3){viewport_width, 0, 0};
    c->vertical = (vec3){0, viewport_height, 0};
    vec3 depth = {0, 0, focal_length};
    c->lower_left_corner = vec3c_sub(vec3c_sub(vec3c_sub(c->origin, vec3_div_k(&(c->horizontal), 2.0)), vec3_div_k(&(c->vertical), 2.0)), depth);
    return c;
}

void camera_free(camera* c){
    free(c);
}

ray camera_get_ray(camera* c, double u, double v){
    vec3 dir = vec3c_sub(vec3c_sum(vec3c_sum(c->lower_left_corner, vec3_mul_k(&c->horizontal, u)), vec3_mul_k(&c->vertical, v)), c->origin);
    return (ray){c->origin, dir};
}



#endif // __CAMERA_H_

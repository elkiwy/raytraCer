#ifndef __CAMERA_H_
#define __CAMERA_H_

//Include standard modules
#include <stdlib.h>

//Include custom modules
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
    vec3 u, v, w;
    double lens_radius;
    double time0, time1;
}camera;


///Initialize the camera object
camera* camera_new(point3 lookfrom, point3 lookat, vec3 vup, double vfov, double aspect_ratio, double aperture, double focus_dist, double time0, double time1){
    ///Alloc memory
    camera* c = malloc(sizeof(camera));


    /*
    ** NOTE:
    ** VFov is the angle of view of the camera,
    ** h is the half height of the image the camera sees with that angle
     */

    //Setup fov and viewport
    double theta = deg2rad(vfov);
    double h = tan(theta/2.0);
    double viewport_height = 2.0*h;
    double viewport_width = aspect_ratio * viewport_height;

    //Setup camera plane
    c->w = vec3c_unit(vec3_sub(&lookfrom, &lookat));
    c->u = vec3c_unit(vec3_cross(&vup, &c->w));
    c->v = vec3_cross(&c->w, &c->u);

    //Setup origin and screens
    c->origin = lookfrom;
    c->horizontal = vec3c_mul_k(vec3_mul_k(&c->u, viewport_width), focus_dist);
    c->vertical = vec3c_mul_k(vec3_mul_k(&c->v, viewport_height), focus_dist);
    c->lower_left_corner = vec3c_sub(vec3c_sub(vec3c_sub(c->origin, vec3_div_k(&(c->horizontal), 2.0)), vec3_div_k(&(c->vertical), 2.0)), vec3_mul_k(&c->w, focus_dist));

    //Setup focus
    c->lens_radius = aperture / 2.0;
    c->time0 = time0;
    c->time1 = time1;
    return c;
}

///Free memory of camera object
void camera_free(camera* c){
    free(c);
}

///Obtain a ray from the camera origin to the specified point, through the lens of the camera to obtain depth of field focus
ray camera_get_ray(camera* c, double s, double t){
    vec3 rd = vec3c_mul_k(random_in_unit_disk(), c->lens_radius);
    vec3 offset = vec3c_sum(vec3_mul_k(&c->u, rd.x), vec3_mul_k(&c->v, rd.y));
    vec3 dir = vec3c_sub(vec3c_sum(vec3c_sum(c->lower_left_corner, vec3_mul_k(&c->horizontal, s)), vec3_mul_k(&c->vertical, t)), c->origin);
    dir = vec3_sub(&dir, &offset);
    vec3 orig = vec3_sum(&c->origin, &offset);
    return (ray){orig, dir, random_double_scaled(c->time0, c->time1)};
}



#endif // __CAMERA_H_

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
camera* camera_new(point3 lookfrom, point3 lookat, vec3 vup, double vfov, double aspect_ratio, double aperture, double focus_dist, double time0, double time1);

///Free memory of camera object
void camera_free(camera* c);

///Obtain a ray from the camera origin to the specified point, through the lens of the camera to obtain depth of field focus
ray camera_get_ray(camera* c, double s, double t);



#endif // __CAMERA_H_

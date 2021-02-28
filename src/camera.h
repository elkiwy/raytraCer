#ifndef __CAMERA_H_
#define __CAMERA_H_

#include <OpenCL/cl.h>

#include "utils.h"




#ifndef PI
#define PI 3.1415926535897932385
#endif


typedef struct camera {
    float origin[3];
    float lower_left_corner[3];
    float horizontal[3];
    float vertical[3];

    cl_float3 w, u, v;
    float lens_radius;
} camera;



camera* make_camera(float3 from, float3 at, float3 vup, float vfov, float aspect_ratio, float aperture, float focus_dist);
void free_camera(camera* cam);

void init_camera(camera* c, cl_float3 from, cl_float3 at, cl_float3 vup, float vfov, float aspect_ratio, float aperture, float focus_dist) ;
void get_ray(const camera* cam, const double s, const double t, cl_float16* ray) ;


#endif // __CAMERA_H_

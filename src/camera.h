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




void init_camera(camera* c, cl_float3 from, cl_float3 at, cl_float3 vup, float vfov, float aspect_ratio, float aperture, float focus_dist) {
    float theta = vfov * PI / 180.0;
    float h = tan(theta/2);
    float viewport_height = 2.0*h;
    float viewport_width = aspect_ratio * viewport_height;

    c->w = float3_normalize((cl_float3){{from.s[0]-at.s[0], from.s[1]-at.s[1], from.s[2]-at.s[2]}});
    c->u = float3_normalize(float3_cross(vup, c->w));
    c->v = float3_cross(c->w, c->u);

    float focal_length = 1.0;
    c->origin[0] = from.s[0];
    c->origin[1] = from.s[1];
    c->origin[2] = from.s[2];
    c->horizontal[0] = focus_dist * viewport_width * c->u.s[0];
    c->horizontal[1] = focus_dist * viewport_width * c->u.s[1];
    c->horizontal[2] = focus_dist * viewport_width * c->u.s[2];
    c->vertical[0] = focus_dist * viewport_height * c->v.s[0];
    c->vertical[1] = focus_dist * viewport_height * c->v.s[1];
    c->vertical[2] = focus_dist * viewport_height * c->v.s[2];
    c->lower_left_corner[0] = c->origin[0] - c->horizontal[0]/2.0 - c->vertical[0]/2.0 - focus_dist*c->w.s[0];
    c->lower_left_corner[1] = c->origin[1] - c->horizontal[1]/2.0 - c->vertical[1]/2.0 - focus_dist*c->w.s[1];
    c->lower_left_corner[2] = c->origin[2] - c->horizontal[2]/2.0 - c->vertical[2]/2.0 - focus_dist*c->w.s[2];

    c->lens_radius = aperture / 2.0;
}


void get_ray(const camera* cam, const double s, const double t, cl_float16* ray) {
    float offset[3];
    cl_float3 rd = random_in_unit_disk();
    rd.s[0] *= cam->lens_radius;
    rd.s[1] *= cam->lens_radius;
    rd.s[2] *= cam->lens_radius;

    offset[0] = cam->u.s[0]*rd.s[0] + cam->v.s[0]*rd.s[1];
    offset[1] = cam->u.s[1]*rd.s[0] + cam->v.s[1]*rd.s[1];
    offset[2] = cam->u.s[2]*rd.s[0] + cam->v.s[2]*rd.s[1];


    /* 0, 1, 2 */
    ray->s[0] = cam->origin[0] + offset[0];
    ray->s[1] = cam->origin[1] + offset[1];
    ray->s[2] = cam->origin[2] + offset[2];
    /* 3, 4, 5 */
    ray->s[3] = cam->lower_left_corner[0] + s*cam->horizontal[0] + t*cam->vertical[0] - cam->origin[0] - offset[0];
    ray->s[4] = cam->lower_left_corner[1] + s*cam->horizontal[1] + t*cam->vertical[1] - cam->origin[1] - offset[1];
    ray->s[5] = cam->lower_left_corner[2] + s*cam->horizontal[2] + t*cam->vertical[2] - cam->origin[2] - offset[2];
    /* 6, 7, 8, 9 */
    ray->s[6] = 0; //dead ray
    ray->s[7] = 0;
    ray->s[8] = 0;
    ray->s[9] = 0;
    /* 10, 11, 12 */ // RG
    ray->s[10] = 1;
    ray->s[11] = 1;
    ray->s[12] = 1;
    /* 13, 14, 15 */
    ray->s[13] = 0;
    ray->s[14] = 0;
    ray->s[15] = 0;
}


#endif // __CAMERA_H_

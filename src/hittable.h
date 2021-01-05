#ifndef __HITTABLE_H_
#define __HITTABLE_H_


#include <stdio.h>
#include <string.h>

#include "vec3.h"
#include "ray.h"
#include "aabb.h"
#include "material.h"
#include "texture.h"


//Forward declare hittable_list
struct hittable_list;




/**
 * Hit record
 * */

//Hit record object to save data about ray hitting hittable object
typedef struct hit_record{
    point3 p;
    vec3 normal;
    struct material* mat;
    double t;
    double u, v;
    int front_face;
}hit_record;






/**
 * Hittables
 * */

///Hittable possible types
typedef enum {
    HITTABLE_BVH_NODE,
    HITTABLE_SPHERE,
    HITTABLE_MOVING_SPHERE,
    HITTABLE_XYRECT,
    HITTABLE_XZRECT,
    HITTABLE_YZRECT,
    HITTABLE_BOX,
    HITTABLE_TRANSLATE,
    HITTABLE_ROTATEY,
    HITTABLE_CONSTANT_MEDIUM,
} hittable_type;

///Generic abstract hittable object
typedef struct hittable{
    hittable_type t;
    void* obj;
} hittable;



/**
 * Hittable normal instances
 * */

///Sphere object
typedef struct { point3 center; double r; struct material* mat; } sphere;
///Axis aligned rect on Z axis
typedef struct { struct material* mat; double x0, x1, y0, y1, k; } xy_rect;
///Axis aligned rect on Y axis
typedef struct { struct material* mat; double x0, x1, z0, z1, k; } xz_rect;
///Axis aligned rect on X axis
typedef struct { struct material* mat; double y0, y1, z0, z1, k; } yz_rect;
///3D Box
typedef struct { point3 box_min; point3 box_max; struct hittable_list* sides; } box;



/**
 * Hittable moving instances
 * */

///Sphere with two centers in 2 times
typedef struct{point3 center0, center1; double time0, time1; double r; struct material* mat;}moving_sphere;



/**
 * Hittable Bounding Volume Hierarchy node
 * */

///Sphere with two centers in 2 times
typedef struct{hittable* left; hittable* right; aabb box;}bvh_node;



/**
 * Hittable Objects wrappers
 * */

///Hittable object wrapper to give constant volume to an object for smoke effects
typedef struct {hittable* boundary; struct material* phase_function; double neg_inv_density;}constant_medium;
///Hittable object Wrapper to apply translation
typedef struct {hittable* obj; vec3 offset;}translate;
///Hittable object Wrapper to apply rotation
//TODO implement all other rotation in a single object
typedef struct {hittable* obj; double sin_theta; double cos_theta; int hasbox; aabb bbox;}rotate_y;





/**
 *
 * Functions
 *
 * */

//Initializers
hittable* hittable_sphere_new(struct hittable_list* world, point3 center, double r, struct material* mat);
hittable* hittable_xy_rect_new(struct hittable_list* world, double x0, double x1, double y0, double y1, double k, struct material* mat);
hittable* hittable_xz_rect_new(struct hittable_list* world, double x0, double x1, double z0, double z1, double k, struct material* mat);
hittable* hittable_yz_rect_new(struct hittable_list* world, double y0, double y1, double z0, double z1, double k, struct material* mat);
hittable* hittable_box_new(struct hittable_list* world, point3 min, point3 max, struct material* mat);
hittable* hittable_moving_sphere_new(struct hittable_list* world, point3 center0, point3 center1, double t0, double t1, double r, struct material* mat);
hittable* bvh_node_init(struct hittable_list* world, struct hittable_list* objects, double t0, double t1);
hittable* hittable_translate_init(struct hittable_list* world, hittable* obj, vec3 offset);
hittable* hittable_rotate_y_init(struct hittable_list* world, hittable* obj, double angle);
hittable* hittable_constant_medium_init(struct hittable_list* world, hittable* b, double d, texture* a);
hittable* hittable_constant_medium_init_c(struct hittable_list* world, hittable* b, double d, color c);

//Deinit
void hittable_free(hittable* o);

//Features
int hittable_hit(hittable* o, ray* r, double tmin, double tmax, hit_record* rec);
int hittable_bounding_box(hittable* o, double t0, double t1, aabb* output_box);


#endif // __HITTABLE_H_

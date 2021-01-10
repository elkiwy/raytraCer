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
    HITTABLE_RECT,
    HITTABLE_BOX,
    HITTABLE_TRANSLATE,
    HITTABLE_ROTATE,
    HITTABLE_CONSTANT_MEDIUM,
    HITTABLE_FLIP_FACE,
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
///Axis aligned rect
typedef enum {XY, XZ, YZ} rect_axis;
typedef struct { struct material* mat; double x0, x1, y0, y1, z0, z1, k; rect_axis axis; } rect;
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
typedef enum {X, Y, Z} rotation_axis;
typedef struct {hittable* obj; double sin_theta; double cos_theta; int hasbox; aabb bbox; rotation_axis axis;}rotate;
///Hittable object Wrapper to flip normals
typedef struct {hittable* obj;}flip_face;



/**
 *
 * Functions
 *
 * */

//Initializers
hittable* hittable_sphere_new(struct hittable_list* world, point3 center, double r, struct material* mat);
hittable* hittable_rect_new(struct hittable_list* world, double x0, double x1, double y0, double y1, double z0, double z1, double k, rect_axis axis, struct material* mat);
hittable* hittable_box_new(struct hittable_list* world, point3 min, point3 max, struct material* mat);
hittable* hittable_moving_sphere_new(struct hittable_list* world, point3 center0, point3 center1, double t0, double t1, double r, struct material* mat);
hittable* bvh_node_init(struct hittable_list* world, struct hittable_list* objects, double t0, double t1);
hittable* hittable_translate_init(struct hittable_list* world, hittable* obj, vec3 offset);
hittable* hittable_rotate_init(struct hittable_list* world, hittable* obj, double angle, rotation_axis axis);
hittable* hittable_constant_medium_init(struct hittable_list* world, hittable* b, double d, texture* a);
hittable* hittable_constant_medium_init_c(struct hittable_list* world, hittable* b, double d, color c);
hittable* hittable_flip_face_init(struct hittable_list* world, hittable* b);

//Deinit
void hittable_free(hittable* o);

//Features
int hittable_hit(hittable* o, ray* r, double tmin, double tmax, hit_record* rec);
int hittable_bounding_box(hittable* o, double t0, double t1, aabb* output_box);

double hittable_pdf_value(hittable* o, point3 orig, vec3 v);
vec3 hittable_random(hittable* o, point3 orig);

#endif // __HITTABLE_H_

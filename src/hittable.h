#ifndef __HITTABLE_H_
#define __HITTABLE_H_


#include "vec3.h"
#include "ray.h"

struct material;

/**
 *
 * Structures
 *
 * */

/** Hit Records */
//Hit record object to save data about ray hitting hittable object
typedef struct{
    point3 p;
    vec3 normal;
    struct material* mat;
    double t;
    int front_face;
}hit_record;


/** Hittables */
//Generic hittable object
typedef enum {HITTABLE_SPHERE} hittable_type;
typedef struct{hittable_type t; void* obj;}hittable;

//Specific hittable sphere object (always stored inside hittable)
typedef struct{
    point3 center;
    double r;
    struct material* mat;
}sphere;


/** Hittable list (world) */
//Container of all the hittable objects
typedef struct{
    hittable** objs;
    int index;
    int size;
}hittable_list;






/**
 *
 * Hit record
 *
 * */

///Set the front_face and normal property of the hit record object based on the ray and normal
inline static void hit_record_set_facenormal(hit_record* rec, ray* r, vec3* outward_normal){
    rec->front_face = vec3_dot(&(r->dir), outward_normal) < 0;
    rec->normal = rec->front_face ? vec3_copy(outward_normal) : vec3_flip(outward_normal);
}

///Print the hit record object
void hit_record_print(hit_record* h){
    printf("Rec: normal: [%.2f, %.2f, %.2f] p: [%.2f, %.2f, %.2f] t: %.2f, front_face: %d\n", h->normal.x, h->normal.y, h->normal.z, h->p.x, h->p.y, h->p.z, h->t, h->front_face );
}





/**
 *
 * Hittable Objects
 *
 * */


///Create a new sphere object
hittable* hittable_sphere_new(hittable_list* world, double x, double y, double z, double r, struct material* mat){
    //Create the sphere
    sphere* s = malloc(sizeof(sphere));
    s->center.x = x;
    s->center.y = y;
    s->center.z = z;
    s->r = r;
    s->mat = mat;

    //Create the hittable generic
    hittable* o = malloc(sizeof(hittable));
    o->obj = s;
    o->t = HITTABLE_SPHERE;

    //Add it to our world
    world->objs[world->index] = o;
    world->index++;
    return o;
}


///Free hittable object
void hittable_free(hittable* o){
    if(o->t == HITTABLE_SPHERE){
        sphere* s = o->obj;
        free(s);
        free(o);
    }else{printf("!! Not implemented hittable_free for type %d\n", o->t);}
}


///Calculate sphere/ray hit and return a hit_record and a success state (0/1)
int sphere_hit(sphere* s, ray* r, double tmin, double tmax, hit_record* rec){
    //Find the discriminant from the formula of intersection between a 3d line and a sphere
    vec3 oc = vec3_sub(&(r->orig), &(s->center));
    double a = vec3_length_squared(&(r->dir));
    double half_b = vec3_dot(&oc, &(r->dir));
    double c = vec3_length_squared(&oc) - s->r * s->r;
    double discriminant = half_b*half_b - a*c;

    //If discriminant is < 0, then there is no intersection...
    if (discriminant < 0){return 0;}
    //...else the sphere and line intersect, so we have to check if the intersection are between tmin and tmax

    //Find nearest root that lies in the acceptable range.
    double sqrtd = sqrt(discriminant);
    double root = (-half_b - sqrtd) / a;
    if (root < tmin || root > tmax){
        //Outside bounds, check the other root
        root = (-half_b + sqrtd) / a;
        if (root < tmin || root > tmax){return 0;}
    }

    //Here we should have valid root, so we fill the data of the hit record
    rec->t = root;
    rec->p = ray_at(r, rec->t);
    vec3 outward_normal = vec3c_div_k(vec3_sub(&rec->p, &s->center), s->r);
    hit_record_set_facenormal(rec, r, &outward_normal);
    rec->mat = s->mat;
    return 1;
}


///Generic function to test for ray hit
int hittable_hit(hittable* o, ray* r, double tmin, double tmax, hit_record* rec){
    if (o->t == HITTABLE_SPHERE){
        return sphere_hit((sphere*)o->obj, r, tmin, tmax, rec);
    }else{printf("!! Not implemented hittable_hit for type %d\n", o->t);return 0;}
}



#endif // __HITTABLE_H_

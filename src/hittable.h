#ifndef __HITTABLE_H_
#define __HITTABLE_H_


#include "vec3.h"
#include "ray.h"

/**
 *
 * Hittable
 *
 * */

typedef struct{
    point3 p;
    vec3 normal;
    double t;
    int front_face;
}hit_record;

inline static void hit_record_set_facenormal(hit_record* rec, ray* r, vec3* outward_normal){
    rec->front_face = (vec3_dot(&(r->dir), outward_normal) < 0);
    rec->normal = rec->front_face ? vec3_copy(outward_normal) : vec3_flip(outward_normal);
}

void hit_record_print(hit_record* h){
    printf("Rec: normal: [%f, %f, %f] p: [%f, %f, %f] t: %f, front_face: %d\n", h->normal.x, h->normal.y, h->normal.z, h->p.x, h->p.y, h->p.z, h->t, h->front_face );
}



typedef enum {HITTABLE_SPHERE} hittable_type;
typedef struct{
    hittable_type t;
    void* obj;
}hittable;

typedef struct{
    point3 center;
    double r;
}sphere;

typedef struct{
    hittable** objs;
    int index;
    int size;
}hittable_list;


hittable* hittable_sphere_new(hittable_list* world, double x, double y, double z, double r){
    //Create the sphere
    sphere* s = malloc(sizeof(sphere));
    s->center.x = x;
    s->center.y = y;
    s->center.z = z;
    s->r = r;

    //Create the hittable generic
    hittable* o = malloc(sizeof(hittable));
    o->obj = s;
    o->t = HITTABLE_SPHERE;

    //Add it to our world
    world->objs[world->index] = o;
    world->index++;

    return o;
}

void hittable_free(hittable* o){
    if(o->t == HITTABLE_SPHERE){
        sphere* s = o->obj;
        free(s);
        free(o);
    }else{printf("!! Not implemented hittable_free for type %d\n", o->t);}
}

///Calculate sphere/ray hit and return a hit_record and a success state (0/1)
int sphere_hit(sphere* s, ray* r, double tmin, double tmax, hit_record* rec){
    vec3 oc = vec3_sub(&(r->orig), &(s->center));
    double a = vec3_length_squared(&(r->dir));
    double half_b = vec3_dot(&oc, &(r->dir));
    double c = vec3_length_squared(&oc) - s->r * s->r;

    double discriminant = half_b*half_b - a*c;
    if (discriminant < 0){return 0;}
    double sqrtd = sqrt(discriminant);

    //Find nearest root that lies in the acceptable range.
    double root = (-half_b - sqrtd) / a;
    if (root < tmin || root > tmax){
        //Outside bounds, check the other root
        root = (-half_b + sqrtd) / a;
        if (root < tmin || root > tmax){return 0;}
    }

    //Here we should have valid root
    rec->t = root;
    rec->p = ray_at(r, rec->t);
    vec3 outward_normal = vec3c_div_k(vec3_sub(&rec->p, &s->center), s->r);
    hit_record_set_facenormal(rec, r, &outward_normal);

    rec->normal = vec3c_div_k(vec3_sub(&(rec->p), &(s->center)), s->r);
    return 1;
}

///Generic function to test for ray hit
int hittable_hit(hittable* o, ray* r, double tmin, double tmax, hit_record* rec){
    if (o->t == HITTABLE_SPHERE){
        return sphere_hit((sphere*)o->obj, r, tmin, tmax, rec);
    }else{printf("!! Not implemented hittable_hit for type %d\n", o->t);return 0;}
}



#endif // __HITTABLE_H_

#ifndef __HITTABLE_H_
#define __HITTABLE_H_


#include <stdio.h>
#include <strings.h>

#include "vec3.h"
#include "ray.h"
#include "aabb.h"

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
    double u, v;
    int front_face;
}hit_record;


/** Hittables */
//Generic hittable object
typedef enum {HITTABLE_BVH_NODE, HITTABLE_SPHERE, HITTABLE_MOVING_SPHERE} hittable_type;
typedef struct{hittable_type t; void* obj;}hittable;

//Specific hittable sphere object (always stored inside hittable)
typedef struct{
    point3 center;
    double r;
    struct material* mat;
}sphere;

//Specific hittable sphere object (always stored inside hittable)
typedef struct{
    point3 center0, center1;
    double time0, time1;
    double r;
    struct material* mat;
}moving_sphere;


typedef struct{
    hittable* left;
    hittable* right;
    aabb box;
}bvh_node;


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
 * Hittable Objects allocators
 *
 * */


///Create a new sphere object
hittable* hittable_sphere_new(hittable_list* world, point3 center, double r, struct material* mat){
    //Create the sphere
    sphere* s = malloc(sizeof(sphere));
    s->center = center;
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

///Create a new moving sphere object
hittable* hittable_moving_sphere_new(hittable_list* world, point3 center0, point3 center1, double t0, double t1, double r, struct material* mat){
    //Create the sphere
    moving_sphere* s = malloc(sizeof(moving_sphere));
    s->center0 = center0;
    s->center1 = center1;
    s->time0 = t0;
    s->time1 = t1;
    s->r = r;
    s->mat = mat;

    //Create the hittable generic
    hittable* o = malloc(sizeof(hittable));
    o->obj = s;
    o->t = HITTABLE_MOVING_SPHERE;

    //Add it to our world
    world->objs[world->index] = o;
    world->index++;
    return o;
}

point3 moving_sphere_center(moving_sphere* s, double t){
    double k = ((t - s->time0) / (s->time1 - s->time0));
    return vec3c_sum(vec3c_mul_k(vec3_sub(&s->center1, &s->center0), k), s->center0);
}














int hittable_bounding_box(hittable* a, double t0, double t1, aabb* output_box);
inline static int box_compare(hittable* a, hittable* b, int axis){
    aabb box_a, box_b;
    if(!hittable_bounding_box(a, 0, 0, &box_a) || !hittable_bounding_box(b, 0, 0, &box_b)){
        printf("Failed to get bounding box in bvh_node comparator.\n"); fflush(stdout);
        return 0;
    }

    if (axis==0){       return box_a.minimum.x < box_b.minimum.x;
    }else if (axis==1){ return box_a.minimum.y < box_b.minimum.y;
    }else{              return box_a.minimum.z < box_b.minimum.z;}
}

int box_x_compare(const void* a, const void* b){return box_compare(*((hittable**)a), *((hittable**)b), 0);}
int box_y_compare(const void* a, const void* b){return box_compare(*((hittable**)a), *((hittable**)b), 1);}
int box_z_compare(const void* a, const void* b){return box_compare(*((hittable**)a), *((hittable**)b), 2);}




hittable* bvh_node_constructor(hittable** src_objects, int src_objects_count, int start, int end, double t0, double t1){
    hittable* h = malloc(sizeof(hittable));
    h->t = HITTABLE_BVH_NODE;

    bvh_node* node = malloc(sizeof(bvh_node));

    hittable** objects = malloc(sizeof(hittable*)*src_objects_count);
    memcpy(objects, src_objects, sizeof(hittable*)*src_objects_count);

    //printf("TEST: %p\n", objects[110]->obj);


    int raxis = random_double();

    int(*comparator)(const void*, const void*);
    if (raxis < 0.33){comparator = &box_x_compare;
    }else if (raxis < 0.66){comparator = &box_y_compare;
    }else{comparator = &box_z_compare;}

    int object_span = end-start;

    if(object_span == 1){
        node->left = objects[start];
        node->right = node->left;

    }else if (object_span == 2){
        if(comparator(&objects[start], &objects[start+1])){
            node->left = objects[start];
            node->right = objects[start+1];
        }else{
            node->left = objects[start+1];
            node->right = objects[start];
        }

    }else{
        qsort(objects, src_objects_count, sizeof(objects[0]), comparator);

        int mid = start + object_span/2;
        node->left  = bvh_node_constructor(src_objects, src_objects_count, start, mid, t0, t1);
        node->right = bvh_node_constructor(src_objects, src_objects_count, mid, end, t0, t1);
    }

    aabb box_left, box_right;
    if(!hittable_bounding_box(node->left, t0, t1, &box_left) || !hittable_bounding_box(node->right, t0, t1, &box_right)){
        printf("Failed to get bounding box in bvh_node constructor.\n"); fflush(stdout);
        return NULL;
    }
    node->box = surrounding_box(box_left, box_right);

    h->obj = node;
    return h;
}


hittable* bvh_node_init(hittable_list* world, hittable_list* objects, double t0, double t1){
    hittable* o = bvh_node_constructor(objects->objs, objects->index, 0, objects->index, t0, t1);

    world->objs[world->index] = o;
    world->index++;

    return o;
}

























/**
 *
 * Hittable Objects deallocators
 *
 * */

///Free hittable object
void hittable_free(hittable* o){
    if(o->t == HITTABLE_SPHERE){
        sphere* s = o->obj;
        free(s);
        free(o);
    }else if(o->t == HITTABLE_MOVING_SPHERE){
        moving_sphere* s = o->obj;
        free(s);
        free(o);
    }else if(o->t == HITTABLE_BVH_NODE){
        bvh_node* n = o->obj;
        free(n);
        free(o);
    }else{printf("!! Not implemented hittable_free for type %d\n", o->t);}
}




/**
 *
 * UV mapping
 *
 * */

void sphere_get_uv(point3* p, double* u, double* v){
    double theta = acos(-p->y);
    double phi = atan2(-p->z, p->x) + PI;
    *u = phi / (2*PI);
    *v = theta / PI;
}




/**
 *
 * Hittable Objects hit testing
 *
 * */
int hittable_hit(hittable* o, ray* r, double tmin, double tmax, hit_record* rec);

int line_sphere_intersection(vec3 oc, vec3 dir, double rad, double tmin, double tmax, double* root){
    double a = vec3_length_squared(&dir);
    double half_b = vec3_dot(&oc, &dir);
    double c = vec3_length_squared(&oc) - rad * rad;
    double discriminant = half_b*half_b - a*c;

    //If discriminant is < 0, then there is no intersection...
    if (discriminant < 0){return 0;}
    //...else the sphere and line intersect, so we have to check if the intersection are between tmin and tmax

    //Find nearest root that lies in the acceptable range.
    double sqrtd = sqrt(discriminant);
    double tmp_root = (-half_b - sqrtd) / a;
    if (tmp_root < tmin || tmp_root > tmax){
        //Outside bounds, check the other root
        tmp_root = (-half_b + sqrtd) / a;
        if (tmp_root < tmin || tmp_root > tmax){
            return 0;
        }
    }

    *root = tmp_root;
    return 1;
}

///Calculate sphere/ray hit and return a hit_record and a success state (0/1)
int sphere_hit(sphere* s, ray* r, double tmin, double tmax, hit_record* rec){
    //Find the discriminant from the formula of intersection between a 3d line and a sphere
    vec3 oc = vec3_sub(&(r->orig), &(s->center));

    double root;
    if (line_sphere_intersection(oc, r->dir, s->r, tmin, tmax, &root)) {
        //Here we should have valid root, so we fill the data of the hit record
        rec->t = root;
        rec->p = ray_at(r, rec->t);
        vec3 outward_normal = vec3c_div_k(vec3_sub(&rec->p, &s->center), s->r);
        hit_record_set_facenormal(rec, r, &outward_normal);
        rec->mat = s->mat;
        sphere_get_uv(&outward_normal, &rec->u, &rec->v);
        return 1;
    }else{
        //Intersection not found
        return 0;
    }
}



///Calculate sphere/ray hit and return a hit_record and a success state (0/1)
int moving_sphere_hit(moving_sphere* s, ray* r, double tmin, double tmax, hit_record* rec){
    //Find the discriminant from the formula of intersection between a 3d line and a sphere
    point3 center_in_time = moving_sphere_center(s, r->time);
    vec3 oc = vec3c_sub(r->orig, center_in_time);

    double root;
    if (line_sphere_intersection(oc, r->dir, s->r, tmin, tmax, &root)) {
        //Here we should have valid root, so we fill the data of the hit record
        rec->t = root;
        rec->p = ray_at(r, rec->t);
        vec3 outward_normal = vec3c_div_k(vec3_sub(&rec->p, &center_in_time), s->r);
        hit_record_set_facenormal(rec, r, &outward_normal);
        rec->mat = s->mat;
        sphere_get_uv(&outward_normal, &rec->u, &rec->v);
        return 1;
    }else{
        //Intersection not found
        return 0;
    }
}



int bvh_node_hit(bvh_node* node, ray* r, double tmin, double tmax, hit_record* rec){
    if(aabb_hit(&node->box, r, tmin, tmax) == 0){return 0;}
    int hit_left = hittable_hit(node->left, r, tmin, tmax, rec);
    int hit_right = hittable_hit(node->right, r, tmin, hit_left ? rec->t : tmax, rec);
    return hit_left || hit_right;
}




///Generic function to test for ray hit
int hittable_hit(hittable* o, ray* r, double tmin, double tmax, hit_record* rec){
    if (o->t == HITTABLE_SPHERE){
        return sphere_hit((sphere*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_MOVING_SPHERE){
        return moving_sphere_hit((moving_sphere*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_BVH_NODE){
        return bvh_node_hit((bvh_node*)o->obj, r, tmin, tmax, rec);
    }else{printf("!! Not implemented hittable_hit for type %d\n", o->t);return 0;}
}





/**
 *
 * Hittable Objects Bounding box
 *
 * */


int sphere_bounding_box(sphere* s, double t0, double t1, aabb* output_box){
    output_box->minimum = vec3c_sub(s->center, (vec3){s->r,s->r,s->r});
    output_box->maximum = vec3c_sum(s->center, (vec3){s->r,s->r,s->r});
    return 1;
}


int moving_sphere_bounding_box(moving_sphere* s, double t0, double t1, aabb* output_box){
    vec3 area = {s->r, s->r, s->r};
    aabb box0 = {vec3c_sub(moving_sphere_center(s, t0), area), vec3c_sum(moving_sphere_center(s, t0), area)};
    aabb box1 = {vec3c_sub(moving_sphere_center(s, t1), area), vec3c_sum(moving_sphere_center(s, t1), area)};
    aabb res = surrounding_box(box0, box1);
    output_box->minimum = res.minimum;
    output_box->maximum = res.maximum;
    return 1;
}



int bvh_node_bounding_box(bvh_node* node, double t0, double t1, aabb* output_box){
    *output_box = node->box;
    return 1;
}

///Generic function to test for ray hit
int hittable_bounding_box(hittable* o, double t0, double t1, aabb* output_box){
    //printf("Bounding box of %p (parent: %p)\n", o->obj, o);fflush(stdout);
    if (o->t == HITTABLE_SPHERE){
        //printf("sphere\n");fflush(stdout);
        return sphere_bounding_box((sphere*)o->obj, t0, t1, output_box);
    }else if (o->t == HITTABLE_MOVING_SPHERE){
        //printf("moving\n");fflush(stdout);
        return moving_sphere_bounding_box((moving_sphere*)o->obj, t0, t1, output_box);
    }else if (o->t == HITTABLE_BVH_NODE){
        //printf("BVH\n");fflush(stdout);
        return bvh_node_bounding_box((bvh_node*)o->obj, t0, t1, output_box);
    }else{printf("!! Not implemented hittable_hit for type %d\n", o->t);return 0;}
}



#endif // __HITTABLE_H_

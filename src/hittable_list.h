#ifndef __HITTABLE_LIST_H_
#define __HITTABLE_LIST_H_

#include <stdlib.h>

#include "vec3.h"
#include "ray.h"
#include "aabb.h"
#include "hittable.h"



/**
 * Hittable list
 * */


//Container of all the hittable objects
typedef struct hittable_list{
    struct hittable** objs;
    int index;
    int size;
    int gc;
}hittable_list;


//Initializers
hittable_list* hittable_list_new_no_gc(int max_size);
hittable_list* hittable_list_new(int max_size);

//Deinits and Garbagecollection
void hittable_list_free(hittable_list* l);
static hittable_list* hittable_lists_allocated[1024];
static int hittable_lists_index = 0;
void hittable_lists_free_all();

//Features
int hittable_list_index(hittable_list* l);
void* hittable_list_objects(hittable_list* l);
void hittable_list_add(hittable_list* l, hittable* o);
int hittable_list_hit(hittable_list* l, ray* r, double tmin, double tmax, hit_record* rec);
int hittable_list_bounding_box(hittable_list* l, double t0, double t1, aabb* output_box);


#endif // __HITTABLE_LIST_H_

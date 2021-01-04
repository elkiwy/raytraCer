#ifndef __HITTABLE_LIST_H_
#define __HITTABLE_LIST_H_

#include <stdlib.h>

#include "vec3.h"
#include "ray.h"
#include "aabb.h"
#include "hittable.h"



/**
 *
 * Hittable list
 *
 * */


/** Hittable list (world) */
//Container of all the hittable objects
typedef struct hittable_list{
    struct hittable** objs;
    int index;
    int size;
}hittable_list;



///Initialize the list object
hittable_list* hittable_list_new(int max_size){
    hittable_list* l = malloc(sizeof(hittable_list));
    l->objs = malloc(sizeof(hittable)*max_size);
    l->index = 0;
    l->size = max_size;
    return l;
}

//Free the list object
void hittable_list_free(hittable_list* l){
    for (int i=0;i<l->index;i++){hittable_free(l->objs[i]);}
    free(l->objs);
    free(l);
}


void hittable_list_add(hittable_list* l, hittable* o){
    if (l==NULL){return;}
    l->objs[l->index] = o;
    l->index++;
}

int hittable_list_index(hittable_list* l){return l->index;}
void* hittable_list_objects(hittable_list* l){return (void*)l->objs;}



///Test a ray to see if it hits any of the objects in the list, and return the hit_record nearest the origin of the camera
int hittable_list_hit(hittable_list* l, ray* r, double tmin, double tmax, hit_record* rec){
    hit_record temp_rec;
    int hit_anything = 0;
    double closest_so_far = tmax;

    //Test all the object between tmin and the current closest hit so far
    for (int i=0; i<l->index; ++i){
        hittable* o = l->objs[i];
        if (hittable_hit(o, r, tmin, closest_so_far, &temp_rec)){
            //Save the hit
            hit_anything = 1;
            closest_so_far = temp_rec.t;

            //Copy current record
            *rec = temp_rec;
        }
    }
    return hit_anything;
}


int hittable_list_bounding_box(hittable_list* l, double t0, double t1, aabb* output_box){
    if (l->index == 0){return 0;}

    aabb temp_box;
    int first_box = 1;

    for (int i=0;i<l->index;++i){
        hittable* o = l->objs[i];
        //Check if everything is alright
        if (hittable_bounding_box(o, t0, t1, &temp_box) == 0){return 0;}
        *output_box = temp_box;
        first_box = 0;
    }

    return 1;
}


#endif // __HITTABLE_LIST_H_

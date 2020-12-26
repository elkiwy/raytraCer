#ifndef __HITTABLE_LIST_H_
#define __HITTABLE_LIST_H_

#include <stdlib.h>

#include "vec3.h"
#include "ray.h"
#include "hittable.h"

/**
 *
 * Hittable list
 *
 * */


hittable_list* hittable_list_new(int max_size){
    hittable_list* l = malloc(sizeof(hittable_list));
    l->objs = malloc(sizeof(hittable)*max_size);
    l->index = 0;
    l->size = max_size;
    return l;
}

void hittable_list_free(hittable_list* l){
    for (int i=0;i<l->index;i++){hittable_free(l->objs[i]);}
    free(l->objs);
    free(l);
}

int hittable_list_hit(hittable_list* l, ray* r, double tmin, double tmax, hit_record* rec){
    hit_record temp_rec;
    int hit_anything = 0;
    double closest_so_far = tmax;

    for (int i=0;i<l->index;++i){
        hittable* o = l->objs[i];
        if (hittable_hit(o, r, tmin, closest_so_far, &temp_rec)){
            hit_anything = 1;
            closest_so_far = temp_rec.t;

            //Copy current record
            *rec = temp_rec;
        }
    }

    return hit_anything;
}


#endif // __HITTABLE_LIST_H_

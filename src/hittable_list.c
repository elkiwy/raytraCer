#include "hittable_list.h"

/**
 *
 * Deinits and garbage collection
 *
 * */


//Free the list object
void hittable_list_free(hittable_list* l){
    //If is marked as GarbageCollected then free all his hittable objects
    if(l->gc){for (int i=0;i<l->index;i++){hittable_free(l->objs[i]);}}
    //Free the list itself
    free(l->objs);
    free(l);
}


///Free all the allocated lists
void hittable_lists_free_all(){
    printf("Freeing %d hittable_lists...\n", hittable_lists_index);fflush(stdout);
    for(int i=0;i<hittable_lists_index;i++){hittable_list_free(hittable_lists_allocated[i]);}
    printf("Done!\n");fflush(stdout);
}





/**
 *
 * Pretty prints
 *
 * */


///Print a single list with all his objects
void hittable_list_print(hittable_list* h){
    printf("\nList %p:\n", h);fflush(stdout);
    for(int i=0;i<h->index;i++){
        printf(" -> %p (%d)\n", h->objs[i], h->objs[i]->t);fflush(stdout);
    }
}


///Prints all the currently allocated lists
void hittable_lists_print_all(){
    for(int i=0;i<hittable_lists_index;i++){
        hittable_list_print(hittable_lists_allocated[i]);
    }
}






/**
 *
 * Initializers
 *
 * */


///Initialize a list object
hittable_list* hittable_list_new(int max_size){
    hittable_list* l = malloc(sizeof(hittable_list));
    l->objs = malloc(sizeof(hittable)*max_size);
    l->index = 0;
    l->size = max_size;
    l->gc = 1;
    hittable_lists_allocated[hittable_lists_index] = l;
    hittable_lists_index++;
    return l;
}


///Initialize the list object without automatic garbage collection of his objects
hittable_list* hittable_list_new_no_gc(int max_size){
    hittable_list* l = hittable_list_new(max_size);
    l->gc = 0;
    return l;
}







/**
 *
 * Features
 *
 * */

///Add a new hittable to a list
void hittable_list_add(hittable_list* l, hittable* o){
    if (l==NULL){return;}
    l->objs[l->index] = o;
    l->index++;
}


///Retrieve the index of a list
int hittable_list_index(hittable_list* l){
    return l->index;
}


///Retrieve the object list
void* hittable_list_objects(hittable_list* l){
    return (void*)l->objs;
}


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


///Test if there is some object whose aabb intersect with t0 and t1
int hittable_list_bounding_box(hittable_list* l, double t0, double t1, aabb* output_box){
    if (l->index == 0){return 0;}

    aabb temp_box;
    int first_box = 1;
    for (int i=0;i<l->index;++i){
        hittable* o = l->objs[i];
        //Check if everything is alright
        if (hittable_bounding_box(o, t0, t1, &temp_box) == 0){return 0;}
        *output_box = first_box ? temp_box : surrounding_box(*output_box, temp_box);
        first_box = 0;
    }

    return 1;
}

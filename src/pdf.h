#ifndef __PDF_H_
#define __PDF_H_

#include "vec3.h"
#include "onb.h"
#include "hittable.h"
#include "hittable_list.h"


//Generic struct
typedef enum {PDF_COSINE, PDF_HITTABLE, PDF_HITTABLE_LIST, PDF_MIXTURE} pdf_type;
typedef struct pdf{pdf_type type; void* obj;}pdf;

//Instances
typedef struct{onb uvw;}pdf_cosine;
typedef struct{point3 orig;hittable* obj;}pdf_hittable;
typedef struct{point3 orig;hittable_list* obj;}pdf_hittable_list;
typedef struct{pdf* a;pdf* b;double ratio;}pdf_mixture;


//Inits
pdf* pdf_generic_init(pdf_type t, void* o);
pdf* pdf_cosine_init(vec3* w);
pdf* pdf_hittable_init(hittable* obj, point3 origin);
pdf* pdf_hittable_list_init(hittable_list* obj, point3 origin);
void pdf_hittable_list_init_stack(pdf* wrapper, pdf_hittable_list* instance, hittable_list* obj, point3 origin);
pdf* pdf_mixture_init(pdf* a, pdf* b, double ratio);
void pdf_mixture_init_stack(pdf* wrapper, pdf_mixture* instance, pdf* a, pdf* b, double ratio);

//Free
void pdf_free(pdf* p);

//Freatures
double pdf_value(pdf* p, vec3* direction);
vec3 pdf_generate(pdf* p);

#endif // __PDF_H_

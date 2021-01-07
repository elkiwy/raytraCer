#ifndef __PDF_H_
#define __PDF_H_

#include "vec3.h"
#include "onb.h"
#include "hittable.h"


typedef enum {PDF_COSINE, PDF_HITTABLE} pdf_type;
typedef struct pdf{
    pdf_type type;
    void* obj;
}pdf;



typedef struct{onb uvw;}pdf_cosine;
typedef struct{point3 orig;hittable* obj;}pdf_hittable;




pdf* pdf_generic_init(pdf_type t, void* o){
    pdf* p = malloc(sizeof(pdf));
    p->type = t;
    p->obj = o;
    return p;
}


pdf* pdf_cosine_init(vec3* w){
    pdf_cosine* c = malloc(sizeof(pdf_cosine));
    onb_init_from_w(&c->uvw, w);
    return pdf_generic_init(PDF_COSINE, c);
}

pdf* pdf_hittable_init(hittable* obj, point3 origin){
    pdf_hittable* h = malloc(sizeof(pdf_hittable));
    h->obj = obj;
    h->orig = origin;
    return pdf_generic_init(PDF_HITTABLE, h);
}






void pdf_free(pdf* p){
    free(p->obj);
    free(p);
}








double pdf_cosine_value(pdf_cosine* p, vec3* dir){
    double cosine = vec3c_dot(vec3_unit(dir), p->uvw.axis[2]);
    return (cosine <= 0) ? 0 : cosine / PI;
}


double pdf_hittable_value(pdf_hittable* p, vec3* dir){
    return hittable_pdf_value(p->obj, p->orig, *dir);
}

double pdf_value(pdf* p, vec3* direction){
    if(p->type == PDF_COSINE){
        return pdf_cosine_value((pdf_cosine*)p->obj, direction);
    }else if(p->type == PDF_HITTABLE){
        return pdf_hittable_value((pdf_hittable*)p->obj, direction);
    }else{
        printf("Not implemented pdf_value");
    }
}







vec3 pdf_cosine_generate(pdf_cosine* p){
    return onb_localv(&p->uvw, random_cosine_direction());
}

vec3 pdf_hittable_generate(pdf_hittable* p){
    return hittable_random(p->obj, p->orig);
}

vec3 pdf_generate(pdf* p){
    if(p->type == PDF_COSINE){
        return pdf_cosine_generate((pdf_cosine*)p->obj);
    }else if(p->type == PDF_HITTABLE){
        return pdf_hittable_generate((pdf_hittable*)p->obj);
    }else{
        printf("Not implemented pdf_generate");

    }
}



#endif // __PDF_H_

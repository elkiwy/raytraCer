#include "pdf.h"



/**
 *
 * Initialize
 *
 * */

///Initialize generic PDF struct
pdf* pdf_generic_init(pdf_type t, void* o){
    pdf* p = malloc(sizeof(pdf));
    p->type = t;
    p->obj = o;
    return p;
}


///Initialize a cosine PDF
pdf* pdf_cosine_init(vec3* w){
    pdf_cosine* c = malloc(sizeof(pdf_cosine));
    onb_init_from_w(&c->uvw, w);
    return pdf_generic_init(PDF_COSINE, c);
}


///Initialize a PDF based on a hittable object
pdf* pdf_hittable_init(hittable* obj, point3 origin){
    pdf_hittable* h = malloc(sizeof(pdf_hittable));
    h->obj = obj;
    h->orig = origin;
    return pdf_generic_init(PDF_HITTABLE, h);
}


///Initialize a PDF based on a hittable list
pdf* pdf_hittable_list_init(hittable_list* obj, point3 origin){
    pdf_hittable_list* h = malloc(sizeof(pdf_hittable_list));
    h->obj = obj;
    h->orig = origin;
    return pdf_generic_init(PDF_HITTABLE_LIST, h);}
///Initialize a PDF based on a hittable list with an existing pdf struct
void pdf_hittable_list_init_stack(pdf* wrapper, pdf_hittable_list* instance, hittable_list* obj, point3 origin){
    instance->obj = obj;
    instance->orig = origin;
    wrapper->obj = instance;
    wrapper->type = PDF_HITTABLE_LIST;
}


///Initialize a mix of two PDFs
pdf* pdf_mixture_init(pdf* a, pdf* b, double ratio){
    pdf_mixture* h = malloc(sizeof(pdf_mixture));
    h->a = a; h->b = b; h->ratio = ratio;
    return pdf_generic_init(PDF_MIXTURE, h);}
///Initialize a mix of two PDFs with an existing pdf struct
void pdf_mixture_init_stack(pdf* wrapper, pdf_mixture* instance, pdf* a, pdf* b, double ratio){
    instance->a = a;
    instance->b = b;
    instance->ratio = ratio;
    wrapper->type = PDF_MIXTURE;
    wrapper->obj = instance;
}





/**
 *
 * Free
 *
 * */

///Free a pdf object
void pdf_free(pdf* p){
    free(p->obj);
    free(p);
}








/**
 *
 * PDF values
 *
 * */

//Forward declaration for recursion
double pdf_value(pdf* p, vec3* direction);


///Value of cosine PDF
double pdf_cosine_value(pdf_cosine* p, vec3* dir){
    double cosine = vec3c_dot(vec3_unit(dir), p->uvw.axis[2]);
    return (cosine <= 0) ? 0 : cosine / PI;
}


///Value of hittable PDF
double pdf_hittable_value(pdf_hittable* p, vec3* dir){
    return hittable_pdf_value(p->obj, p->orig, *dir);
}


///Value of hittable list PDF
double pdf_hittable_list_value(pdf_hittable_list* p, vec3* dir){
    return hittable_list_pdf_value(p->obj, &p->orig, dir);
}


///Value of mixture PDF
double pdf_mixture_value(pdf_mixture* p, vec3* dir){
    return p->ratio * pdf_value(p->a, dir) + (1-p->ratio) * pdf_value(p->b, dir);
}


///Generic PDF value function call
double pdf_value(pdf* p, vec3* direction){
    if(p->type == PDF_COSINE){
        return pdf_cosine_value((pdf_cosine*)p->obj, direction);
    }else if(p->type == PDF_HITTABLE){
        return pdf_hittable_value((pdf_hittable*)p->obj, direction);
    }else if(p->type == PDF_HITTABLE_LIST){
        return pdf_hittable_list_value((pdf_hittable_list*)p->obj, direction);
    }else if(p->type == PDF_MIXTURE){
        return pdf_mixture_value((pdf_mixture*)p->obj, direction);
    }else{
        printf("Not implemented pdf_value");
    }
}








/**
 *
 * Random point from PDF
 *
 * */

//Forward declaration for recursion
vec3 pdf_generate(pdf* p);


///Generate random point of cosine PDF
vec3 pdf_cosine_generate(pdf_cosine* p){
    return onb_localv(&p->uvw, random_cosine_direction());
}


///Generate random point of hittable PDF
vec3 pdf_hittable_generate(pdf_hittable* p){
    return hittable_random(p->obj, p->orig);
}


///Generate random point of hittable list PDF
vec3 pdf_hittable_list_generate(pdf_hittable_list* p){
    return hittable_list_random(p->obj, p->orig);
}


///Generate random point of mixture PDF
vec3 pdf_mixture_generate(pdf_mixture* p){
    if (random_double() < p->ratio){
        return pdf_generate(p->a);
    }else{
        return pdf_generate(p->b);
    }
}


///Generic PDF generate function call
vec3 pdf_generate(pdf* p){
    if(p->type == PDF_COSINE){
        return pdf_cosine_generate((pdf_cosine*)p->obj);
    }else if(p->type == PDF_HITTABLE){
        return pdf_hittable_generate((pdf_hittable*)p->obj);
    }else if(p->type == PDF_HITTABLE_LIST){
        return pdf_hittable_list_generate((pdf_hittable_list*)p->obj);
    }else if(p->type == PDF_MIXTURE){
        return pdf_mixture_generate((pdf_mixture*)p->obj);
    }else{
        printf("Not implemented pdf_generate");
        return (vec3){0,0,0};
    }
}

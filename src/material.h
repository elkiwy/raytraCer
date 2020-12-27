#ifndef __MATERIAL_H_
#define __MATERIAL_H_

#include <stdio.h>
#include <stdlib.h>

#include "vec3.h"
#include "ray.h"
#include "hittable.h"




//Generic material struct
typedef enum{MATERIAL_LAMBERTIAN, MATERIAL_METAL, MATERIAL_DIELECTRIC} material_type;
typedef struct{
    material_type type;
    void* mat;
}material;

//Materials
typedef struct{color albedo;}material_lambertian;
typedef struct{color albedo; double fuzz;}material_metal;
typedef struct{double ir;}material_dielectric;



material* material_lambertian_new(color a){
    material* m = malloc(sizeof(material));
    m->type = MATERIAL_LAMBERTIAN;
    material_lambertian* l = malloc(sizeof(material_lambertian));
    l->albedo = a;
    m->mat = l;
    return m;
}

material* material_metal_new(color a, double fuzz){
    material* m = malloc(sizeof(material));
    m->type = MATERIAL_METAL;
    material_metal* metal = malloc(sizeof(material_metal));
    metal->albedo = a;
    metal->fuzz = fuzz;
    m->mat = metal;
    return m;
}


material* material_dielectric_new(double ir){
    material* m = malloc(sizeof(material));
    m->type = MATERIAL_DIELECTRIC;
    material_dielectric* d = malloc(sizeof(material_dielectric));
    d->ir = ir;
    m->mat = d;
    return m;
}




int material_lambertian_scatter(material_lambertian* mat, ray*r, hit_record* rec, color* attenuation, ray* scattered){
    vec3 scatter_direction = vec3c_sum(rec->normal, vec3_random_unit_vector_in_unit_sphere());
    //Check for near zero scattered directions
    if (vec3_is_near_zero(&scatter_direction)){scatter_direction = rec->normal;}
    *scattered = (ray){rec->p, scatter_direction};
    *attenuation = mat->albedo;
    return 1;
}


int material_metal_scatter(material_metal* mat, ray*r, hit_record* rec, color* attenuation, ray* scattered){
    vec3 reflected = vec3c_reflect(vec3_unit(&r->dir), rec->normal);
    vec3 scattered_direction = vec3c_sum(reflected, vec3c_mul_k(vec3_random_in_unit_sphere(), mat->fuzz));
    *scattered = (ray){rec->p, scattered_direction};
    *attenuation = mat->albedo;
    double a = vec3_dot(&scattered->dir, &rec->normal);
    return a > 0;
}



double reflectance(double cosine, double ref_idx){
    double r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0*r0;
    return r0 + (1-r0)*pow((1-cosine),5);
}

int material_dielectric_scatter(material_dielectric* mat, ray*r, hit_record* rec, color* attenuation, ray* scattered){
    *attenuation = (color){1,1,1};
    double refraction_ratio = rec->front_face ? (1.0/mat->ir) : mat->ir;
    vec3 unit_dir = vec3_unit(&r->dir);

    //Check if the ray can be refracted or it needs to be reflected
    double cos_theta = fmin(vec3c_dot(vec3_flip(&unit_dir), rec->normal), 1.0);
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    int cannot_refract = refraction_ratio * sin_theta > 1.0;
    vec3 direction;
    if(cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double()){
        direction = vec3_reflect(&unit_dir, &rec->normal);
    }else{
        direction = vec3c_refract(unit_dir, rec->normal, refraction_ratio);
    }

    *scattered = (ray){rec->p, direction};
    return 1;
}





int material_scatter(material* mat, ray* r, hit_record* rec, color* attenuation, ray* scattered){
    if(mat->type == MATERIAL_LAMBERTIAN){
        return material_lambertian_scatter((material_lambertian*)mat->mat, r, rec, attenuation, scattered);
    }else if(mat->type == MATERIAL_DIELECTRIC){
        return material_dielectric_scatter((material_dielectric*)mat->mat, r, rec, attenuation, scattered);
    }else if(mat->type == MATERIAL_METAL){
        return material_metal_scatter((material_metal*)mat->mat, r, rec, attenuation, scattered);
    }else{printf("material_scatter nor implemented for material type %d\n", mat->type); return 0;}
}









#endif // __MATERIAL_H_

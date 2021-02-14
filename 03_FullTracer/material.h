#ifndef __MATERIAL_H_
#define __MATERIAL_H_

#include <stdio.h>
#include <stdlib.h>

#include "vec3.h"
#include "ray.h"
#include "hittable.h"
#include "texture.h"
#include "onb.h"

//Forward declaration of hit record
struct hit_record;
struct pdf;

///Holder structure for data abour scattering rays
typedef struct scatter_record{
    ray specular_ray;
    int is_specular;
    color attenuation;
    struct pdf* pdf_ptr;
}scatter_record;


//Generic material struct
typedef enum{MATERIAL_LAMBERTIAN, MATERIAL_METAL, MATERIAL_DIELECTRIC, MATERIAL_LIGHT, MATERIAL_ISOTROPIC} material_type;
typedef struct material {material_type type; void* mat;}material;

//Materials
typedef struct{texture* albedo;}material_lambertian;
typedef struct{color albedo; double fuzz;}material_metal;
typedef struct{double ir;}material_dielectric;
typedef struct{texture* emit;}material_light;
typedef struct{texture* albedo;}material_isotropic;


//Init
material* material_lambertian_new(texture* t);
material* material_lambertian_new_from_color(color c);
material* material_metal_new(color a, double fuzz);
material* material_dielectric_new(double ir);
material* material_light_new(texture* t);
material* material_light_new_from_color(color c);
material* material_isotropic_new(texture* t);
material* material_isotropic_new_from_color(color c);

//Deinit
void material_free(material* m);
static material* materials_allocated[1024];
static int materials_index = 0;
void materials_free_all();

//Features
color material_emitted(material* mat, struct hit_record* rec, double u, double v, point3 p);
int material_scatter(material* mat, ray* r, struct hit_record* rec, scatter_record* srec);
double material_scattering_pdf(material* mat, ray* r, struct hit_record* rec, ray* scattered);

#endif // __MATERIAL_H_

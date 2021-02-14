#include "material.h"




/**
 * Deinits and Garbage collection
 * */

///Free a material
void material_free(material* m){
    free(m->mat);
    free(m);
}


///Free all the allocated materials
void materials_free_all(){
    printf("Freeing %d materials...\n", materials_index);fflush(stdout);
    for(int i=0;i<materials_index;i++){material_free(materials_allocated[i]);}
    printf("Done!\n");fflush(stdout);
}





/**
 * Materials Initializations
 */


///Init a generic abstract material
material* material_generic_init(material_type t){
    material* m = malloc(sizeof(material));
    m->type = t;

    //Add to allocated array
    materials_allocated[materials_index] = m;
    materials_index++;
    return m;
}


///Init Lambertian
material* material_lambertian_new(texture* t){
    material* m = material_generic_init(MATERIAL_LAMBERTIAN);
    material_lambertian* l = malloc(sizeof(material_lambertian));
    l->albedo = t;
    m->mat = l;
    return m;}
///Init Lambertian with a color
material* material_lambertian_new_from_color(color c){
    texture* t = texture_solid_color_init(c);
    return material_lambertian_new(t);
}


///Init metal
material* material_metal_new(color a, double fuzz){
    material* m = material_generic_init(MATERIAL_METAL);
    material_metal* metal = malloc(sizeof(material_metal));
    metal->albedo = a;
    metal->fuzz = fuzz;
    m->mat = metal;
    return m;
}


///Init dielectric
material* material_dielectric_new(double ir){
    material* m = material_generic_init(MATERIAL_DIELECTRIC);
    material_dielectric* d = malloc(sizeof(material_dielectric));
    d->ir = ir;
    m->mat = d;
    return m;
}


///Init light
material* material_light_new(texture* t){
    material* m = material_generic_init(MATERIAL_LIGHT);
    material_light* l = malloc(sizeof(material_light));
    l->emit = t;
    m->mat = l;
    return m;}
///Init light with color
material* material_light_new_from_color(color c){
    texture* t = texture_solid_color_init(c);
    return material_light_new(t);
}


///Init isotropic
material* material_isotropic_new(texture* t){
    material* m = material_generic_init(MATERIAL_ISOTROPIC);
    material_isotropic* l = malloc(sizeof(material_isotropic));
    l->albedo = t;
    m->mat = l;
    return m;}
///Init isotropic with color
material* material_isotropic_new_from_color(color c){
    texture* t = texture_solid_color_init(c);
    return material_isotropic_new(t);
}








/**
 * Material Scattering
 */


///Scatter a ray that hits a lambertian material
int material_lambertian_scatter(material_lambertian* mat, ray*r, hit_record* rec, color* attenuation, ray* scattered){
    //Scatter the ray following the normal of the hit, and slightly randomize it
    vec3 scatter_direction = vec3c_sum(rec->normal, vec3_random_unit_vector_in_unit_sphere());
    //Check for near zero scattered directions
    if (vec3_is_near_zero(&scatter_direction)){scatter_direction = rec->normal;}
    *scattered = (ray){rec->p, scatter_direction, r->time};
    *attenuation = texture_value(mat->albedo, rec->u, rec->v, &rec->p);
    return 1;
}


///Scatter a ray that hits a metal material
int material_metal_scatter(material_metal* mat, ray*r, hit_record* rec, color* attenuation, ray* scattered){
    //Scatter the ray reflecting it using the normal as bisector
    vec3 reflected = vec3c_reflect(vec3_unit(&r->dir), rec->normal);
    //Randomize it with fuzz
    vec3 scattered_direction = vec3c_sum(reflected, vec3c_mul_k(vec3_random_in_unit_sphere(), mat->fuzz));
    *scattered = (ray){rec->p, scattered_direction, r->time};
    *attenuation = mat->albedo;
    double a = vec3_dot(&scattered->dir, &rec->normal);
    return a > 0;
}


///Reflectance formula
double reflectance(double cosine, double ref_idx){double r0 = (1-ref_idx) / (1+ref_idx); r0 = r0*r0; return r0 + (1-r0)*pow((1-cosine),5);}
///Scatter a ray that hits a dielectric material
int material_dielectric_scatter(material_dielectric* mat, ray*r, hit_record* rec, color* attenuation, ray* scattered){
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
        direction = vec3c_refract(unit_dir, rec->normal, refraction_ratio);}
    *scattered = (ray){rec->p, direction, r->time};
    *attenuation = (color){1,1,1};
    return 1;
}


///Scatter a light
int material_light_scatter(){
    return 0;
}


///Scatter a ray that hits a isotropic material
int material_isotropic_scatter(material_isotropic* mat, ray*r, hit_record* rec, color* attenuation, ray* scattered){
    *scattered = (ray){rec->p, vec3_random_in_unit_sphere(), r->time};
    *attenuation = texture_value(mat->albedo, rec->u, rec->v, &rec->p);
    return 1;
}


///Generic scatter material
int material_scatter(material* mat, ray* r, hit_record* rec, color* attenuation, ray* scattered){
    if(mat->type == MATERIAL_LAMBERTIAN){
        return material_lambertian_scatter((material_lambertian*)mat->mat, r, rec, attenuation, scattered);
    }else if(mat->type == MATERIAL_DIELECTRIC){
        return material_dielectric_scatter((material_dielectric*)mat->mat, r, rec, attenuation, scattered);
    }else if(mat->type == MATERIAL_METAL){
        return material_metal_scatter((material_metal*)mat->mat, r, rec, attenuation, scattered);
    }else if(mat->type == MATERIAL_LIGHT){
        return material_light_scatter();
    }else if(mat->type == MATERIAL_ISOTROPIC){
        return material_isotropic_scatter((material_isotropic*)mat->mat, r, rec, attenuation, scattered);
    }else{printf("material_scatter nor implemented for material type %d\n", mat->type); return 0;}
}







/**
 * Material Emission
 */

///Light emission
color material_light_emitted(material_light* light, double u, double v, point3 p){
    return texture_value(light->emit, u, v, &p);
}


///Material emission
color material_emitted(material* mat, double u, double v, point3 p){
    //Return the emitted color
    if(mat->type == MATERIAL_LIGHT){
        return material_light_emitted((material_light*)mat->mat, u, v, p);
    //Else is not emmitting
    }else{return (color){0,0,0};}
}

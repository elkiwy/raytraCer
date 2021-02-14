#include "ray.h"

#include "hittable.h"
#include "hittable_list.h"
#include "pdf.h"

///Move the origin point through the ray for T distance
inline point3 ray_at(ray* r, double t){
    vec3 a = vec3_mul_k(&r->dir, t);
    return vec3_sum(&r->orig, &a);
}

///Print a ray object
void ray_print(char* s, ray* r){
    printf("%s: orig: [%.2f %.2f %.2f] dir: [%.2f %.2f %.2f]\n", s, r->orig.x, r->orig.y, r->orig.z, r->dir.x, r->dir.y, r->dir.z);
}






/**
 *
 * Main
 *
 * */

///Prints a color in hex RGB
void color_printHex(char* s, color* c){
    printf("%s: %02x %02x %02x\n", s, (int)(256*clamp(c->x, 0.0, 0.999)), (int)(256*clamp(c->y, 0.0, 0.999)), (int)(256*clamp(c->z, 0.0, 0.999)));
}

///Write a pixel color to the output file, the color is averaged between the samples_per_pixel
void write_color(unsigned char* pixels, color* pixel_color, int samples_per_pixel, int index){
    double r = pixel_color->x; double g = pixel_color->y; double b = pixel_color->z;
    double scale = 1.0 / (double)samples_per_pixel;

    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);


    pixels[index] = (unsigned char)(256*clamp(r, 0.0, 0.999));;
    pixels[index+1] = (unsigned char)(256*clamp(g, 0.0, 0.999));;
    pixels[index+2] = (unsigned char)(256*clamp(b, 0.0, 0.999));;
}

///Cast a ray into the world and retrieve the color of that ray
color ray_color(ray* r, color* background, hittable_list* world, hittable_list* lights, int recursion_depth){
    //Base case of recursion
    if (recursion_depth <= 0){return (color){0,0,0};}

    //If the ray doesn't hits anything, take the color of the background
    hit_record rec;
    if (hittable_list_hit(world, r, 0.001, HUGE_VAL, &rec) == 0){return *background;}

    //If the thing I hit is emitting something, take the emitted color
    scatter_record srec;
    color emitted = material_emitted(rec.mat, &rec, rec.u, rec.v, rec.p);

    //If the hit wasn't scattered take the emitted color (could be also just a black if it wasn't emitting anything)
    if (material_scatter(rec.mat, r, &rec, &srec) == 0){return emitted;}
    if (srec.is_specular){return vec3c_mul(srec.attenuation, ray_color(&srec.specular_ray, background, world, lights, recursion_depth-1));}

    //Get the lights PDF to generate the scattered ray
    pdf light_pdf; pdf_hittable_list light_instance;
    pdf_hittable_list_init_stack(&light_pdf, &light_instance, lights, rec.p);

    //Mix the perfect light PDF with a random cosine one
    pdf mixture_pdf; pdf_mixture mixture_instance;
    pdf_mixture_init_stack(&mixture_pdf, &mixture_instance, &light_pdf, srec.pdf_ptr, 0.1);

    //Generate the scattered ray
    ray scattered = {rec.p, pdf_generate(&mixture_pdf), r->time};
    double pdf_val = pdf_value(&mixture_pdf, &scattered.dir);
    if(srec.pdf_ptr != 0){pdf_free(srec.pdf_ptr);}

    //else Recursively take the color of the scattered ray and multiply its color by the attenuation color
    color scattered_pdf_color = vec3_mul_k(&srec.attenuation,  material_scattering_pdf(rec.mat, r, &rec, &scattered));
    color scattered_ray_color = ray_color(&scattered, background, world, lights, recursion_depth-1);
    color scattered_ray_color_weighted_by_pdf = vec3_div_k(&scattered_ray_color, pdf_val);
    return vec3c_sum(emitted, vec3_mul(&scattered_pdf_color, &scattered_ray_color_weighted_by_pdf));
}

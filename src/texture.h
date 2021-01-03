#ifndef __TEXTURE_H_
#define __TEXTURE_H_


#include "extern_stb_image.h"

#include "utils.h"
#include "vec3.h"
#include "perlin.h"





typedef enum {TEXTURE_SOLID, TEXTURE_CHECKER, TEXTURE_NOISE, TEXTURE_IMAGE} texture_type;


typedef struct {
    texture_type type;
    void* o;
} texture;

typedef struct{
    color color_value;
}solid_color;

typedef struct{
    texture* odd;
    texture* even;
}checker;

typedef struct{
    perlin* noise;
    double scale;
}texture_noise;

typedef struct{
    unsigned char* data;
    int width, height;
    int bytes_per_scanline;
}texture_image;



/**
 *
 * Initializers
 *
 * */

texture* texture_solid_color_init(color c){
    texture* t = malloc(sizeof(texture));
    t->type = TEXTURE_SOLID;
    solid_color* sc = malloc(sizeof(solid_color));
    sc->color_value = c;
    t->o = sc;
    return t;}
texture* texture_solid_color_init_rgb(double r, double g, double b){
    return texture_solid_color_init((color){r, g, b});
}



texture* texture_checker_init(texture* odd, texture* even){
    texture* t = malloc(sizeof(texture));
    t->type = TEXTURE_CHECKER;
    checker* chk = malloc(sizeof(checker));
    chk->even = even;
    chk->odd = odd;
    t->o = chk;
    return t;}
texture* texture_checker_init_c(color odd, color even){
    texture* odd_t = texture_solid_color_init(odd);
    texture* even_t = texture_solid_color_init(even);
    return texture_checker_init(odd_t, even_t);
}



texture* texture_noise_init_scaled(double sc){
    texture* t = malloc(sizeof(texture));
    t->type = TEXTURE_NOISE;
    texture_noise* n = malloc(sizeof(texture_noise));
    perlin* p = perlin_init();
    n->noise = p;
    n->scale = sc;
    t->o = n;
    return t;}
texture* texture_noise_init(){
    return texture_noise_init_scaled(1.0);
}



#define BYTES_PER_PIXEL 3
texture* texture_image_init(char* filename){
    texture* t = malloc(sizeof(texture));
    t->type = TEXTURE_IMAGE;
    texture_image* i = malloc(sizeof(texture_image));
    int components_per_pixel = BYTES_PER_PIXEL;
    i->data = stbi_load(filename, &i->width, &i->height, &components_per_pixel, components_per_pixel);
    if(!i->data){
        printf("!! could not load texture image file '%s'! \n", filename);fflush(stdout);
        return NULL;
    }

    i->bytes_per_scanline = BYTES_PER_PIXEL * i->width;
    t->o = i;
    return t;
}





/**
 *
 * Deallocators
 *
 * */

void texture_free(texture* t){
    if (t->type == TEXTURE_SOLID){
        free(t->o);
        free(t);
    }else if (t->type == TEXTURE_CHECKER){
        free(t->o);
        free(t);
    }else if (t->type == TEXTURE_NOISE){
        free(t->o);
        free(t);
    }else if (t->type == TEXTURE_IMAGE){
        free(((texture_image*)(t->o))->data);
        free(t->o);
        free(t);
    }else{printf("texture_free not implemented for type %d\n", t->type);fflush(stdout);}
}



/**
 *
 * Color values
 *
 * */

color texture_value(texture* t, double u, double v, vec3* p);

color solid_color_value(solid_color* sc, double u, double v, vec3* p){
    return sc->color_value;
}

color checker_value(checker* sc, double u, double v, vec3* p){
    double sines = sin(10*p->x) * sin(10*p->y) * sin(10*p->z);
    return (sines < 0) ? texture_value(sc->odd, u, v, p) : texture_value(sc->even, u, v, p);
}

color noise_texture_value(texture_noise* tn, double u, double v, vec3* p){
    //double perl = perlin_turb(tn->noise, vec3_mul_k(p, tn->scale));
    double perl = perlin_turb(tn->noise, *p);
    double val = 0.5 * (1.0 + sin(tn->scale*p->z + 10 * perl));
    return vec3c_mul_k((color){1,1,1}, val);
}


color image_texture_value(texture_image* img, double u, double v, vec3* p){
    if(!img->data){return (color){1,0,1};}

    u = clamp(u, 0.0, 1.0);
    v = 1.0 - clamp(v, 0.0, 1.0);

    int i = (int)(u*img->width);
    int j = (int)(v*img->height);

    if (i>=img->width)  i = img->width - 1;
    if (j>=img->height) j = img->height - 1;

    const double color_scale = 1.0 / 255.0;
    unsigned char* pixel = img->data + j*img->bytes_per_scanline + i*BYTES_PER_PIXEL;

    return (color){color_scale*pixel[0], color_scale*pixel[1], color_scale*pixel[2]};
}




color texture_value(texture* t, double u, double v, vec3* p){
    if (t->type == TEXTURE_SOLID){
        return solid_color_value((solid_color*)t->o, u, v, p);
    }else if (t->type == TEXTURE_CHECKER){
        return checker_value((checker*)t->o, u, v, p);
    }else if (t->type == TEXTURE_NOISE){
        return noise_texture_value((texture_noise*)t->o, u, v, p);
    }else if (t->type == TEXTURE_IMAGE){
        return image_texture_value((texture_image*)t->o, u, v, p);
    }else{printf("texture_value not implemented for type %d\n", t->type);fflush(stdout); return (color){0,0,0};}
}






#endif // __TEXTURE_H_

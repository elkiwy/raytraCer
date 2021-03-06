#ifndef __TEXTURE_H_
#define __TEXTURE_H_



#include "utils.h"
#include "vec3.h"
#include "perlin.h"




/**
 * Generic struc
 * */

typedef enum {TEXTURE_SOLID, TEXTURE_CHECKER, TEXTURE_NOISE, TEXTURE_IMAGE} texture_type;
typedef struct {texture_type type; void* o;} texture;


/**
 * Instances
 * */

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

//Init
texture* texture_solid_color_init(color c);
texture* texture_solid_color_init_rgb(double r, double g, double b);
texture* texture_checker_init(texture* odd, texture* even);
texture* texture_checker_init_c(color odd, color even);
texture* texture_noise_init_scaled(double sc);
texture* texture_noise_init();
texture* texture_image_init(char* filename);

//Free
void texture_free(texture* t);
static texture* textures_allocated[1024];
static int textures_index = 0;
void textures_free_all();

//Feature
color texture_value(texture* t, double u, double v, vec3* p);

#endif // __TEXTURE_H_

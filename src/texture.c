#include "texture.h"


#include "extern_stb_image.h"
#define BYTES_PER_PIXEL 3


/**
 *
 * Initializers
 *
 * */

///Init generic object
texture* texture_generic_init(texture_type type, void* obj){
    texture* t = malloc(sizeof(texture));
    t->type = type;
    t->o = obj;
    textures_allocated[textures_index] = t;
    textures_index++;
    return t;
}


///Solid color texture from color
texture* texture_solid_color_init(color c){
    solid_color* sc = malloc(sizeof(solid_color));
    sc->color_value = c;
    return texture_generic_init(TEXTURE_SOLID, sc);}
///Solid color texture from RGB values
texture* texture_solid_color_init_rgb(double r, double g, double b){
    return texture_solid_color_init((color){r, g, b});
}


///Checker texture with two textures
texture* texture_checker_init(texture* odd, texture* even){
    checker* chk = malloc(sizeof(checker));
    chk->even = even;
    chk->odd = odd;
    return texture_generic_init(TEXTURE_CHECKER, chk);}
///Checker texture with two colors
texture* texture_checker_init_c(color odd, color even){
    texture* odd_t = texture_solid_color_init(odd);
    texture* even_t = texture_solid_color_init(even);
    return texture_checker_init(odd_t, even_t);
}


///Noise texture with scale
texture* texture_noise_init_scaled(double sc){
    texture_noise* n = malloc(sizeof(texture_noise));
    perlin* p = perlin_init();
    n->noise = p;
    n->scale = sc;
    return texture_generic_init(TEXTURE_NOISE, n);}
///Noise texture
texture* texture_noise_init(){
    return texture_noise_init_scaled(1.0);
}


///Image Texture init
texture* texture_image_init(char* filename){
    texture_image* i = malloc(sizeof(texture_image));
    int components_per_pixel = BYTES_PER_PIXEL;
    i->data = stbi_load(filename, &i->width, &i->height, &components_per_pixel, components_per_pixel);
    if(!i->data){printf("!! could not load texture image file '%s'! \n", filename);fflush(stdout); return NULL;}
    i->bytes_per_scanline = BYTES_PER_PIXEL * i->width;
    return texture_generic_init(TEXTURE_IMAGE, i);
}





/**
 *
 * Deallocators
 *
 * */

///Free a texture
void texture_free(texture* t){
    if (t->type == TEXTURE_NOISE){perlin_free(((texture_noise*)t->o)->noise);
    }else if (t->type == TEXTURE_IMAGE){free(((texture_image*)(t->o))->data);}
    free(t->o);
    free(t);
}


///Free all the allocated textures
void textures_free_all(){
    printf("Freeing %d textures...\n", textures_index);fflush(stdout);
    for(int i=0;i<textures_index;i++){texture_free(textures_allocated[i]);}
    printf("Done!\n");fflush(stdout);
}







/**
 *
 * Color values
 *
 * */

//Forward declare texture_value for cyclic
color texture_value(texture* t, double u, double v, vec3* p);


///Solid color value
color solid_color_value(solid_color* sc){
    return sc->color_value;
}


///Checker value
color checker_value(checker* sc, double u, double v, vec3* p){
    double sines = sin(10*p->x) * sin(10*p->y) * sin(10*p->z);
    return (sines < 0) ? texture_value(sc->odd, u, v, p) : texture_value(sc->even, u, v, p);
}


///Noise value
color noise_texture_value(texture_noise* tn, vec3* p){
    double perl = perlin_turb(tn->noise, *p);
    double val = 0.5 * (1.0 + sin(tn->scale*p->z + 10 * perl));
    return vec3c_mul_k((color){1,1,1}, val);
}


///Image value
color image_texture_value(texture_image* img, double u, double v){
    //Default color for no image
    if(!img->data){return (color){1,0,1};}

    //Map UV to image
    u = clamp(u, 0.0, 1.0);
    v = 1.0 - clamp(v, 0.0, 1.0);
    int i = (int)(u*img->width);
    int j = (int)(v*img->height);
    if (i>=img->width)  i = img->width - 1;
    if (j>=img->height) j = img->height - 1;

    //Scale color and find the right pixel
    const double color_scale = 1.0 / 255.0;
    unsigned char* pixel = img->data + j*img->bytes_per_scanline + i*BYTES_PER_PIXEL;
    return (color){color_scale*pixel[0], color_scale*pixel[1], color_scale*pixel[2]};
}


///Main texture value function call
color texture_value(texture* t, double u, double v, vec3* p){
    if (t->type == TEXTURE_SOLID){
        return solid_color_value((solid_color*)t->o);
    }else if (t->type == TEXTURE_CHECKER){
        return checker_value((checker*)t->o, u, v, p);
    }else if (t->type == TEXTURE_NOISE){
        return noise_texture_value((texture_noise*)t->o, p);
    }else if (t->type == TEXTURE_IMAGE){
        return image_texture_value((texture_image*)t->o, u, v);
    }else{printf("texture_value not implemented for type %d\n", t->type);fflush(stdout); return (color){0,0,0};}
}

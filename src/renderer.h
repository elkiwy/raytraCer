#ifndef __RENDERER_H_
#define __RENDERER_H_


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/time.h>

#include <OpenCL/cl.h>

//Include OpemMP for multithreading
#include <omp.h>



//STB Image write to output png
#include "external/stb_image_write.h"

#include "utils.h"
#include "camera.h"
#include "objects.h"





typedef struct scene_settings{
    int res_width;
    int res_height;
    int samples_per_pixel;
    int iterations;
}scene_settings;

typedef struct scene{
    cl_float16* objs;
    cl_float16* lights;
    cl_float16* wrapped_objs;
    cl_float16* mats;
    cl_float16* texs;

    int objs_max;
    int lights_max;
    int wrapped_objs_max;
    int mats_max;
    int texs_max;

    int objs_count;
    int lights_count;
}scene;


scene* make_scene(const ObjectList* objects, const ObjectList* lights, const int size);
void free_scene(scene* s);



int get_optimal_chunk_splitting(const int img_w, const int img_h, const int samples_per_pixel, cl_device_id dev);

void render_scene(const scene* s, const camera* cam, const scene_settings* settings, const char* filename);


#endif // __RENDERER_H_

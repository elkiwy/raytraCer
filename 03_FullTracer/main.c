//Include standard libs
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//Include OpemMP for multithreading
#include <omp.h>

//Include every module
#include "vec3.h"
#include "ray.h"
#include "hittable.h"
#include "hittable_list.h"
#include "camera.h"
#include "utils.h"
#include "material.h"
#include "pdf.h"

//STB Image write to output png
#include "external/stb_image_write.h"






/**
 *
 * Scene
 *
 * */

///Cornel box
hittable_list* setupScene(hittable_list* lights_output){
    hittable_list* world = hittable_list_new(1024);

    material* red = material_lambertian_new_from_color((color){0.65, 0.05, 0.05});
    material* whi = material_lambertian_new_from_color((color){0.73, 0.73, 0.73});
    material* gre = material_lambertian_new_from_color((color){0.12, 0.45, 0.15});
    material* lig = material_light_new_from_color((color){15,15,15});

    //Box and light
    hittable* b1 = hittable_rect_new(world, 0,0, 0, 555, 0, 555, 555, YZ, gre);
    hittable* b2 = hittable_rect_new(world, 0,0, 0, 555, 0, 555, 0, YZ, red);
    hittable* b3 = hittable_rect_new(NULL, 213, 343, 0,0, 227, 332, 554, XZ, lig);
    b3 = hittable_flip_face_init(world, b3);
    hittable_list_add(lights_output, b3);
    hittable* b4 = hittable_rect_new(world, 0, 555, 0,0, 0, 555, 0, XZ, whi);
    hittable* b5 = hittable_rect_new(world, 0, 555, 0,0, 0, 555, 555, XZ, whi);
    hittable* b6 = hittable_rect_new(world, 0, 555, 0, 555, 0,0, 555, XY, whi);


    material* aluminium = material_metal_new((color){0.8, 0.85, 0.88}, 0.0);
    hittable* box1 = hittable_box_new(NULL, (point3){0,0,0}, (point3){165,330,165}, aluminium);
    box1 = hittable_rotate_init(NULL, box1, 15, Y);
    box1 = hittable_translate_init(world, box1, (vec3){265, 0, 295});

    material* glass = material_dielectric_new(1.5);
    hittable* s2 = hittable_sphere_new(world, (point3){190, 90, 190}, 90, glass);
    hittable_list_add(lights_output, s2);

    return world;
}



hittable_list* simpleScene(hittable_list* lights_output){
    hittable_list* world = hittable_list_new(1024);
    material* red = material_lambertian_new_from_color((color){0.65, 0.05, 0.05});
    material* lig = material_light_new_from_color((color){15,15,15});

    //Box and light
    hittable* b3 = hittable_rect_new(NULL, 213, 343, 0,0, 227, 332, 554, XZ, lig);
    b3 = hittable_flip_face_init(world, b3);
    hittable_list_add(lights_output, b3);

    hittable* s2 = hittable_sphere_new(world, (point3){190, 190, 190}, 90, red);

    return world;
}



///Main
int main(int argc, char** argv) {
    //Read parameters
    char* outputFilePath = "";
    double ASPECT_RATIO = 1.0;
    int IMAGE_WIDTH = 1024;
    int IMAGE_HEIGHT = (int)(IMAGE_WIDTH / ASPECT_RATIO);
    int samples_per_pixel = 256; //how many ray per pixel
    int max_recursion_depth = 10; //how deep the ray scattering goes
    for (int i = 1; i < argc ; i++) {
        if (argv[i][0] == '-'){
            switch (argv[i][1]) {
            case 'w': IMAGE_WIDTH  = atoi(argv[i+1]); break;
            case 'h': IMAGE_HEIGHT = atoi(argv[i+1]); break;
            case 'o': outputFilePath = argv[i+1]; break;
            case 's': samples_per_pixel = atoi(argv[i+1]); break;
            default:
                fprintf(stderr, "Usage: %s [-w width] [-h height] [-o outputFile] [-s sample_per_pixel]\n", argv[0]);
                exit(EXIT_FAILURE);
            }
        }
    }

    //Define camera
    double vfov = 40;
    vec3 vup = (vec3){0, 1, 0};
    double dist_to_focus = 10;
    double aperture = 0.0;
    point3 lookfrom = (point3){278, 278, -800};
    point3 lookat   = (point3){278, 278, 0};

    //Define objects list
    hittable_list* lights = hittable_list_new_no_gc(64);
    //hittable_list* world = setupScene(lights);
    hittable_list* world = simpleScene(lights);
    color background = {0,0,0};

    //Init camera
    camera* c = camera_new(lookfrom, lookat, vup, vfov, ASPECT_RATIO, aperture, dist_to_focus, 0.0, 1.0);

    //PNG output setup
    const int CHANNELS = 3;
    unsigned char pixels[IMAGE_WIDTH * IMAGE_HEIGHT * CHANNELS];

    //Setup profiling and multithreading
    double time = 0.0;
    double begin = omp_get_wtime();
    int progress = 0;
    const unsigned int THREADS = 4*4;
    int step = floor(IMAGE_HEIGHT/THREADS);
    int chunksDone = 0;

    //Split height into chunks and give them to threads
    #pragma omp parallel for
    for(unsigned int k=0; k<THREADS; ++k){
        double thread_begin = omp_get_wtime();

        //Each threads cycle a step numer of rows
        for(int sj=step-1; sj>=0; --sj){
            int j = sj+(k*step);
            //Cycle pixel inside each row
            for(int i=0; i<IMAGE_WIDTH; ++i){
                //Sample the same pixel with variations on the casted ray to create antialiasing
                color pixel_color = {0,0,0};
                for(int s=0; s<samples_per_pixel; ++s){
                    //Cast ray and sum the color
                    double u = ((double)i+random_double()) / (double)(IMAGE_WIDTH - 1);
                    double v = ((double)j+random_double()) / (double)(IMAGE_HEIGHT - 1);
                    ray r = camera_get_ray(c, u, v);
                    pixel_color = vec3c_sum(pixel_color, ray_color(&r, &background, world, lights, max_recursion_depth));
                }

                //Output the result color into the file
                int index = (i+(j*IMAGE_WIDTH)) * CHANNELS;
                write_color(pixels, &pixel_color, samples_per_pixel, index);
            }
        }

        //Output feedback of chunk completed
        int thread_id = omp_get_thread_num();
        double thread_end = omp_get_wtime();
        chunksDone++;
        printf("Thread %d (chunk: %d-%d) finished in %f, remains %d chunks.\n", thread_id, k*step, (k+1)*step, (double)(thread_end - thread_begin), THREADS-chunksDone);
    }

    //Flip the final image and save it
    stbi_flip_vertically_on_write(1);
    stbi_write_png(outputFilePath, IMAGE_WIDTH, IMAGE_HEIGHT, CHANNELS, pixels, IMAGE_WIDTH * CHANNELS);

    //Output total time elapsed
    double end = omp_get_wtime();
    time = (double)(end - begin);
    printf("Time elpased for rendering %f\n", time);

    //Free memory
    camera_free(c);
    hittable_lists_free_all();
    materials_free_all();
    textures_free_all();

    //Exit
    return 0;
}

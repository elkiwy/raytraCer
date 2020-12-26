//Include standard libs
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//Include every module
#include "vec3.h"
#include "ray.h"
#include "hittable.h"
#include "hittable_list.h"
#include "camera.h"
#include "utils.h"


/**
 *
 * Main
 *
 * */

///Write a pixel color to the output file, the color is averaged between the samples_per_pixel
void write_color(FILE* f, color* pixel_color, int samples_per_pixel){
    double r = pixel_color->x; double g = pixel_color->y; double b = pixel_color->z;
    double scale = 1.0 / (double)samples_per_pixel;
    r *= scale; g *= scale; b *= scale;
    fprintf(f, "%d %d %d\n", (int)(256*clamp(r, 0.0, 0.999)), (int)(256*clamp(g, 0.0, 0.999)), (int)(256*clamp(b, 0.0, 0.999)));
}


///Cast a ray into the world and retrieve the color of that ray
color ray_color(ray* r, hittable_list* world){
    //Define background
    color white = {1,1,1};
    color cyan = {0.5, 0.7, 1.0};

    //If the ray hits anything, take the color of the object
    hit_record rec;
    if (hittable_list_hit(world, r, 0, HUGE_VAL, &rec)){
        return vec3c_mul_k(vec3_sum(&(rec.normal), &white), 0.5);
    }

    //Else draw a pixel of the background white/cyan gradient
    vec3 unit_dir = vec3_unit(&r->dir);
    double t = 0.5 * (unit_dir.y + 1.0);
    return vec3c_sum(vec3_mul_k(&white, (1.0 - t)), vec3_mul_k(&cyan, t));
}



///Main
int main(int argc, char** argv) {
    //Output file
    char* outputFilePath = argv[1];
    char buffer[1024*1024]; // 1 MB buffer
    FILE* f = fopen(outputFilePath, "w+");

    //Define camera
    camera* c = camera_new();
    const double ASPECT_RATIO = 16.0/9.0;
    const int IMAGE_WIDTH = 400;
    const int IMAGE_HEIGHT = (int)(IMAGE_WIDTH / ASPECT_RATIO);
    const int samples_per_pixel = 32;

    //Define objects list
    hittable_list* world = hittable_list_new(32);
    hittable_sphere_new(world, 0,      0, -1, 0.5);
    hittable_sphere_new(world, 0, -100.5, -1, 100);

    //Processing
    int progress = 0;
    fprintf(f, "P3\n%d %d\n255\n", IMAGE_WIDTH, IMAGE_HEIGHT);
    //Cycle Rows of pixels
    for(int j=IMAGE_HEIGHT-1; j>=0; --j){
        //Print progress
        if (progress < (1-((double)j/(double)(IMAGE_HEIGHT-1)))*100){printf("Progress: %d%%\n",progress); progress += 10;}

        //Cycle pixel inside each row
        for(int i=0; i<IMAGE_WIDTH; ++i){
            //Sample the same pixel with variations on the casted ray to create antialiasing
            color pixel_color = {0,0,0};
            for(int s=0; s<samples_per_pixel; ++s){
                double u = ((double)i+random_double()) / (double)(IMAGE_WIDTH - 1);
                double v = ((double)j+random_double()) / (double)(IMAGE_HEIGHT - 1);
                ray r = camera_get_ray(c, u, v);
                pixel_color = vec3c_sum(pixel_color, ray_color(&r, world));
            }

            //Output the result color into the file
            write_color(f, &pixel_color, samples_per_pixel);
        }
    }

    //Complete progress
    printf("Progress: 100%%\n");
    fclose(f);

    //Free memory
    camera_free(c);
    hittable_list_free(world);

    //Exit
    return 0;
}

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
#include "material.h"


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
void write_color(FILE* f, color* pixel_color, int samples_per_pixel){
    double r = pixel_color->x; double g = pixel_color->y; double b = pixel_color->z;
    double scale = 1.0 / (double)samples_per_pixel;
    r = sqrt(scale * r); g = sqrt(scale * g); b = sqrt(scale * b);
    fprintf(f, "%d %d %d\n", (int)(256*clamp(r, 0.0, 0.999)), (int)(256*clamp(g, 0.0, 0.999)), (int)(256*clamp(b, 0.0, 0.999)));
}


///Cast a ray into the world and retrieve the color of that ray
color ray_color(ray* r, hittable_list* world, int recursion_depth){
    //Base case of recursion
    if (recursion_depth <= 0){return (color){0,0,0};}

    //If the ray hits anything, take the color of the object
    hit_record rec;
    if (hittable_list_hit(world, r, 0.001, HUGE_VAL, &rec)){
        //Check if the thing the ray hit scatter the ray somewhere else or absorbs the ray
        ray scattered;
        color attenuation;
        if (material_scatter((material*)rec.mat, r, &rec, &attenuation, &scattered)){
            //Recursively take the color of the scattered ray and multiply its color by the attenuation color
            color scattered_ray_color = ray_color(&scattered, world, recursion_depth-1);
            return vec3_mul(&attenuation, &scattered_ray_color);
        }else{
            //Ray got absorbed by the material
            return (color){0,0,0};
        }
    }

    //Else draw a pixel of the background white/cyan gradient
    color white = {1,1,1};
    color cyan = {0.5, 0.7, 1.0};
    vec3 unit_dir = vec3_unit(&r->dir);
    double t = 0.5 * (unit_dir.y + 1.0);
    return vec3c_sum(vec3_mul_k(&white, (1.0 - t)), vec3_mul_k(&cyan, t));
}





///Setup all the objects in the scene
hittable_list* setup_scene(){
    //Create the world container
    hittable_list* world = hittable_list_new(1024);

    //Add the ground object
    material* material_ground = material_lambertian_new((color){0.5, 0.5, 0.5});
    hittable_sphere_new(world,  0, -1000, 0, 1000,  (struct material*)material_ground);

    //Add small random spheres
    for (int a=-11; a<11; a++){
        for (int b=-11; b<11; b++){
            double choose_mat = random_double();
            point3 center = {a + 0.9*random_double(), 0.2, b + 0.9*random_double()};
            vec3 p = {4, 0.2, 0};
            if ( vec3c_length(vec3_sub(&center, &p)) > 0.9){
                material* sphere_mat;
                if(choose_mat < 0.8){
                    //diffuse
                    color albedo = vec3c_mul(vec3_random() , vec3_random());
                    sphere_mat = material_lambertian_new(albedo);
                    hittable_sphere_new(world, center.x, center.y, center.z, 0.2, (struct material*)sphere_mat);
                }else if(choose_mat < 0.95){
                    //metal
                    color albedo = vec3_random_scaled(0.5, 1);
                    double fuzz = random_double_scaled(0, 0.5);
                    sphere_mat = material_metal_new(albedo, fuzz);
                    hittable_sphere_new(world, center.x, center.y, center.z, 0.2, (struct material*)sphere_mat);
                }else{
                    //glass
                    sphere_mat = material_dielectric_new(1.5);
                    hittable_sphere_new(world, center.x, center.y, center.z, 0.2, (struct material*)sphere_mat);
                }
            }
        }
    }

    //Add 3 big shperes with different materials and return
    material* material1 = material_dielectric_new(1.5);
    hittable_sphere_new(world, 0, 1, 0, 1.0, (struct material*)material1);
    material* material2 = material_lambertian_new((color){0.4,0.2,0.1});
    hittable_sphere_new(world, -4, 1, 0, 1.0, (struct material*)material2);
    material* material3 = material_metal_new((color){0.7, 0.6, 0.5}, 0.0);
    hittable_sphere_new(world,  4, 1, 0, 1.0, (struct material*)material3);
    return world;
}












///Main
int main(int argc, char** argv) {
    /**
     *
     * Setup
     *
     **/
    //Output file
    char* outputFilePath = argv[1];
    char buffer[1024*1024]; // 1 MB buffer
    FILE* outputFile = fopen(outputFilePath, "w+");

    //Define camera
    const double ASPECT_RATIO = 16.0/9.0;
    point3 lookfrom = (point3){13, 2,  3};
    point3 lookat   = (point3){ 0, 0,  0};
    vec3 vup        = (vec3){0, 1, 0};
    double dist_to_focus = 10;
    double aperture = 0.1;
    camera* c = camera_new(lookfrom, lookat, vup, 30, ASPECT_RATIO, aperture, dist_to_focus);

    //Define render target
    const int IMAGE_WIDTH = 1920;
    const int IMAGE_HEIGHT = (int)(IMAGE_WIDTH / ASPECT_RATIO);
    const int samples_per_pixel = 8; //how many ray per pixel
    const int max_recursion_depth = 8; //how deep the ray scattering goes

    //Define objects list
    hittable_list* world = setup_scene();



    /**
     *
     * Rendering
     *
     **/

    //Initialize output PPM
    fprintf(outputFile, "P3\n%d %d\n255\n", IMAGE_WIDTH, IMAGE_HEIGHT);

    //Cycle Rows of pixels
    int progress = 0;
    for(int j=IMAGE_HEIGHT-1; j>=0; --j){
        //Print progress
        if (progress < (1-((double)j/(double)(IMAGE_HEIGHT-1)))*100){printf("Progress: %d%%\n",progress); progress += 10;}

        //Cycle pixel inside each row
        for(int i=0; i<IMAGE_WIDTH; ++i){
            //Sample the same pixel with variations on the casted ray to create antialiasing
            color pixel_color = {0,0,0};
            for(int s=0; s<samples_per_pixel; ++s){
                //Cast ray and sum the color
                double u = ((double)i+random_double()) / (double)(IMAGE_WIDTH - 1);
                double v = ((double)j+random_double()) / (double)(IMAGE_HEIGHT - 1);
                ray r = camera_get_ray(c, u, v);
                pixel_color = vec3c_sum(pixel_color, ray_color(&r, world, max_recursion_depth));
            }

            //Output the result color into the file
            write_color(outputFile, &pixel_color, samples_per_pixel);
        }
    }

    //Complete progress
    printf("Progress: 100%%\n");
    fclose(outputFile);

    //Free memory
    camera_free(c);
    hittable_list_free(world);

    //Exit
    return 0;
}

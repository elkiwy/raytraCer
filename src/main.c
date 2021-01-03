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
color ray_color(ray* r, color* background, hittable_list* world, int recursion_depth){
    //Base case of recursion
    if (recursion_depth <= 0){return (color){0,0,0};}

    //If the ray hits anything, take the color of the object
    hit_record rec;
    if (hittable_list_hit(world, r, 0.001, HUGE_VAL, &rec) == 0){return *background;}

    //Check if the thing the ray hit scatter the ray somewhere else or absorbs the ray
    ray scattered;
    color attenuation;
    color emitted = material_emitted((material*)rec.mat, rec.u, rec.v, rec.p);

    if (material_scatter((material*)rec.mat, r, &rec, &attenuation, &scattered) == 0){return emitted;}

    //Recursively take the color of the scattered ray and multiply its color by the attenuation color
    color scattered_ray_color = ray_color(&scattered, background, world, recursion_depth-1);
    return vec3c_sum(emitted, vec3_mul(&attenuation, &scattered_ray_color));
}




///Setup all the objects in the scene
hittable_list* setup_scene(){
    //Create the world container
    hittable_list* world = hittable_list_new(1024);

    texture* ground_texture1 = texture_solid_color_init_rgb(0.5, 0.5, 0.5);
    texture* ground_texture2 = texture_solid_color_init_rgb(0.2, 0.2, 0.2);
    texture* ground_texture = texture_checker_init(ground_texture1, ground_texture2);
    texture* ball_texture = texture_solid_color_init_rgb(0.4, 0.2, 0.1);

    //Add the ground object
    material* material_ground = material_lambertian_new(ground_texture);
    hittable* g = hittable_sphere_new(world,  (point3){0, -1000, 0}, 1000,  (struct material*)material_ground);

    //Add small random spheres
    //for (int a=-11; a<11; a++){
    //    for (int b=-11; b<11; b++){
    //        double choose_mat = random_double();
    //        point3 center = {a + 0.9*random_double(), 0.2, b + 0.9*random_double()};
    //        vec3 p = {4, 0.2, 0};
    //        if ( vec3c_length(vec3_sub(&center, &p)) > 0.9){
    //            material* sphere_mat;
    //            if(choose_mat < 0.8){
    //                //diffuse
    //                color albedo = vec3c_mul(vec3_random() , vec3_random());
    //                sphere_mat = material_lambertian_new(albedo);


    //                point3 center2 = vec3c_sum(center, (vec3){0, random_double_scaled(0, 0.5), 0});
    //                hittable_moving_sphere_new(world, center, center2, 0.0, 1.0, 0.2, (struct material*)sphere_mat);

    //            }else if(choose_mat < 0.95){
    //                //metal
    //                color albedo = vec3_random_scaled(0.5, 1);
    //                double fuzz = random_double_scaled(0, 0.5);
    //                sphere_mat = material_metal_new(albedo, fuzz);
    //                hittable_sphere_new(world, center, 0.2, (struct material*)sphere_mat);
    //            }else{
    //                //glass
    //                sphere_mat = material_dielectric_new(1.5);
    //                hittable_sphere_new(world, center, 0.2, (struct material*)sphere_mat);
    //            }
    //        }
    //    }
    //}

    //Add 3 big shperes with different materials and return
    material* material1 = material_dielectric_new(1.5);
    hittable* h1 = hittable_sphere_new(world, (point3){0, 1, 0}, 1.0, (struct material*)material1);
    material* material2 = material_lambertian_new(ball_texture);
    hittable* h2 = hittable_sphere_new(world, (point3){-4, 1, 0}, 1.0, (struct material*)material2);
    material* material3 = material_metal_new((color){0.7, 0.6, 0.5}, 0.0);
    hittable* h3 = hittable_sphere_new(world,  (point3){4, 1, 0}, 1.0, (struct material*)material3);

    printf("ground: %p \n", g);
    printf("palle: %p %p %p\n", h1, h2, h3);
    printf("world: \n"); for(int i=0;i<world->index;i++){printf("==> world[%i]: %p (pointer at %p)\n", i, world->objs[i], &(world->objs[i]) );}


    hittable_list* world2 = hittable_list_new(1024);
    bvh_node_init(world2, world, 0, 1);


    printf("world2: \n"); for(int i=0;i<world2->index;i++){printf("==> world2[%i]: %p (pointer at %p)\n", i, world2->objs[i], &(world2->objs[i]) );}

    return world2;
}








hittable_list* two_spheres(){
    hittable_list* world = hittable_list_new(1024);
    texture* texture1 = texture_solid_color_init_rgb(0.5, 0.5, 0.5);
    texture* texture2 = texture_solid_color_init_rgb(0.2, 0.2, 0.2);
    texture* texture = texture_checker_init(texture1, texture2);
    material* mat1 = material_lambertian_new(texture);
    hittable* s1 = hittable_sphere_new(world, (point3){0, -10, 0}, 10, (struct material*)mat1);
    hittable* s2 = hittable_sphere_new(world, (point3){0,  10, 0}, 10, (struct material*)mat1);
    return world;
}



hittable_list* two_perlin_spheres(){
    hittable_list* world = hittable_list_new(1024);

    texture* pertext = texture_noise_init_scaled(4);
    material* mat1 = material_lambertian_new(pertext);
    material* mat2 = material_lambertian_new(pertext);
    hittable* s1 = hittable_sphere_new(world, (point3){0, -1000, 0}, 1000, (struct material*)mat1);
    hittable* s2 = hittable_sphere_new(world, (point3){0,  2, 0}, 2, (struct material*)mat1);
    return world;
}


hittable_list* earth(){
    hittable_list* world = hittable_list_new(1024);

    texture* earthtexture = texture_image_init("images/earthmap.jpg");
    material* earthsurface = material_lambertian_new(earthtexture);
    hittable* globe = hittable_sphere_new(world, (point3){0, 0, 0}, 2, (struct material*)earthsurface);
    return world;
}



hittable_list* simple_light(){
    hittable_list* world = hittable_list_new(1024);

    texture* pertext = texture_noise_init_scaled(4);
    material* mat1 = material_lambertian_new(pertext);
    hittable* s1 = hittable_sphere_new(world, (point3){0, -1000, 0}, 1000, (struct material*)mat1);
    hittable* s2 = hittable_sphere_new(world, (point3){0,  2, 0}, 2, (struct material*)mat1);

    material* difflight = material_light_new_from_color((color){4,4,4});
    hittable* l1 = hittable_xy_rect_new(world, 3, 5, 1, 3, -2, (struct material*)difflight);

    return world;
}



hittable_list* cornell_box(){
    hittable_list* world = hittable_list_new(1024);

    material* red = material_lambertian_new_from_color((color){0.65, 0.05, 0.05});
    material* whi = material_lambertian_new_from_color((color){0.73, 0.73, 0.73});
    material* gre = material_lambertian_new_from_color((color){0.12, 0.45, 0.15});
    material* lig = material_light_new_from_color((color){15,15,15});


    hittable* b1 = hittable_yz_rect_new(world, 0, 555, 0, 555, 555, (struct material*)gre);
    hittable* b2 = hittable_yz_rect_new(world, 0, 555, 0, 555, 0, (struct material*)red);

    hittable* b3 = hittable_xz_rect_new(world, 213, 343, 227, 332, 554, (struct material*)lig);
    hittable* b4 = hittable_xz_rect_new(world, 0, 555, 0, 555, 0, (struct material*)whi);
    hittable* b5 = hittable_xz_rect_new(world, 0, 555, 0, 555, 555, (struct material*)whi);

    hittable* b6 = hittable_xy_rect_new(world, 0, 555, 0, 555, 555, (struct material*)whi);
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
    printf("Argc: %d\n", argc);
    char* outputFilePath = argv[1];
    char buffer[1024*1024]; // 1 MB buffer
    FILE* outputFile = fopen(outputFilePath, "w+");

    //Define render target
    double ASPECT_RATIO = 16.0/9.0;
    //Define camera
    double vfov = 40;
    point3 lookfrom;
    point3 lookat;
    vec3 vup = (vec3){0, 1, 0};
    double dist_to_focus = 10;
    double aperture = 0.0;
    color background = {0,0,0};

    //Define objects list
    hittable_list* world;
    switch(0){
        case 1:
            world = setup_scene();
            background = (color){0.70, 0.80, 1.00};
            lookfrom = (point3){13, 2,  3};
            lookat   = (point3){ 0, 0,  0};
            aperture = 0.1;
            break;

        case 2:
            world = two_spheres();
            background = (color){0.70, 0.80, 1.00};
            lookfrom = (point3){13, 2,  3};
            lookat   = (point3){ 0, 0,  0};
            vfov = 20.0;
            break;

        case 3:
            world = two_perlin_spheres();
            background = (color){0.70, 0.80, 1.00};
            lookfrom = (point3){13, 2,  3};
            lookat   = (point3){ 0, 0,  0};
            vfov = 20.0;
            break;

        case 4:
            world = earth();
            background = (color){0.70, 0.80, 1.00};
            lookfrom = (point3){13, 2,  3};
            lookat   = (point3){ 0, 0,  0};
            vfov = 20.0;
            break;

        case 5:
            world = simple_light();
            background = (color){0.0, 0.0, 0.0};
            lookfrom = (point3){26, 3, 6};
            lookat   = (point3){ 0, 2, 0};
            vfov = 20.0;
            break;

        default:
        case 6:
            world = cornell_box();
            ASPECT_RATIO = 1.0;
            background = (color){0.0, 0.0, 0.0};
            lookfrom = (point3){278, 278, -800};
            lookat   = (point3){278, 278, 0};
            vfov = 40.0;
            break;
    }

    //Init camera
    const int IMAGE_WIDTH = 320;
    const int IMAGE_HEIGHT = (int)(IMAGE_WIDTH / ASPECT_RATIO);
    const int samples_per_pixel = 128; //how many ray per pixel
    const int max_recursion_depth = 8; //how deep the ray scattering goes
    camera* c = camera_new(lookfrom, lookat, vup, vfov, ASPECT_RATIO, aperture, dist_to_focus, 0.0, 1.0);


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
                pixel_color = vec3c_sum(pixel_color, ray_color(&r, &background, world, max_recursion_depth));
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

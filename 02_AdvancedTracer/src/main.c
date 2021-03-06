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

//STB Image write to output png
#include "extern_stb_image_write.h"

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
    r = sqrt(scale * r); g = sqrt(scale * g); b = sqrt(scale * b);
    pixels[index] = (unsigned char)(256*clamp(r, 0.0, 0.999));;
    pixels[index+1] = (unsigned char)(256*clamp(g, 0.0, 0.999));;
    pixels[index+2] = (unsigned char)(256*clamp(b, 0.0, 0.999));;
}


///Cast a ray into the world and retrieve the color of that ray
color ray_color(ray* r, color* background, hittable_list* world, int recursion_depth){
    //Base case of recursion
    if (recursion_depth <= 0){return (color){0,0,0};}

    //If the ray doesn't hits anything, take the color of the background
    hit_record rec;
    if (hittable_list_hit(world, r, 0.001, HUGE_VAL, &rec) == 0){return *background;}
    ray scattered;
    color attenuation;

    //If the thing I hit is emitting something, take the emitted color
    color emitted = material_emitted((material*)rec.mat, rec.u, rec.v, rec.p);

    //If the hit wasn't scattered take the emitted color (could be also just a black if it wasn't emitting anything)
    if (material_scatter((material*)rec.mat, r, &rec, &attenuation, &scattered) == 0){return emitted;}

    //else Recursively take the color of the scattered ray and multiply its color by the attenuation color
    color scattered_ray_color = ray_color(&scattered, background, world, recursion_depth-1);
    return vec3c_sum(emitted, vec3_mul(&attenuation, &scattered_ray_color));
}





/**
 *
 * Scenes
 *
 * */

///3 big balls and many smaller ones
hittable_list* setup_scene(){
    hittable_list* world = hittable_list_new(1024);

    //Checker texture
    texture* ground_texture1 = texture_solid_color_init_rgb(0.5, 0.5, 0.5);
    texture* ground_texture2 = texture_solid_color_init_rgb(0.2, 0.2, 0.2);
    texture* ground_texture = texture_checker_init(ground_texture1, ground_texture2);
    texture* ball_texture = texture_solid_color_init_rgb(0.4, 0.2, 0.1);

    //Add the ground object
    material* material_ground = material_lambertian_new(ground_texture);
    hittable* g = hittable_sphere_new(world,  (point3){0, -1000, 0}, 1000,  (struct material*)material_ground);

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
                    sphere_mat = material_lambertian_new_from_color(albedo);

                    point3 center2 = vec3c_sum(center, (vec3){0, random_double_scaled(0, 0.5), 0});
                    hittable_moving_sphere_new(world, center, center2, 0.0, 1.0, 0.2, (struct material*)sphere_mat);

                }else if(choose_mat < 0.95){
                    //metal
                    color albedo = vec3_random_scaled(0.5, 1);
                    double fuzz = random_double_scaled(0, 0.5);
                    sphere_mat = material_metal_new(albedo, fuzz);
                    hittable_sphere_new(world, center, 0.2, (struct material*)sphere_mat);
                }else{
                    //glass
                    sphere_mat = material_dielectric_new(1.5);
                    hittable_sphere_new(world, center, 0.2, (struct material*)sphere_mat);
                }
            }
        }
    }

    //Add 3 big shperes with different materials and return
    material* material1 = material_dielectric_new(1.5);
    hittable* h1 = hittable_sphere_new(world, (point3){0, 1, 0}, 1.0, (struct material*)material1);
    material* material2 = material_lambertian_new(ball_texture);
    hittable* h2 = hittable_sphere_new(world, (point3){-4, 1, 0}, 1.0, (struct material*)material2);
    material* material3 = material_metal_new((color){0.7, 0.6, 0.5}, 0.0);
    hittable* h3 = hittable_sphere_new(world,  (point3){4, 1, 0}, 1.0, (struct material*)material3);
    return world;
}


///Two big spheres with checker textures
hittable_list* two_spheres(){
    hittable_list* world = hittable_list_new(1024);
    texture* texture1 = texture_solid_color_init_rgb(0.5, 0.5, 0.5);
    texture* texture2 = texture_solid_color_init_rgb(0.2, 0.2, 0.2);
    texture* tex = texture_checker_init(texture1, texture2);
    material* mat1 = material_lambertian_new(tex);
    hittable* s1 = hittable_sphere_new(world, (point3){0, -10, 0}, 10, (struct material*)mat1);
    hittable* s2 = hittable_sphere_new(world, (point3){0,  10, 0}, 10, (struct material*)mat1);
    return world;
}


///One sphere and ground with perlin noise
hittable_list* two_perlin_spheres(){
    hittable_list* world = hittable_list_new(1024);
    texture* pertext = texture_noise_init_scaled(4);
    material* mat1 = material_lambertian_new(pertext);
    material* mat2 = material_lambertian_new(pertext);
    hittable* s1 = hittable_sphere_new(world, (point3){0, -1000, 0}, 1000, (struct material*)mat1);
    hittable* s2 = hittable_sphere_new(world, (point3){0,  2, 0}, 2, (struct material*)mat1);
    return world;
}


///Earth sphere image
hittable_list* earth(){
    hittable_list* world = hittable_list_new(1024);
    texture* earthtexture = texture_image_init("images/earthmap.jpg");
    material* earthsurface = material_lambertian_new(earthtexture);
    hittable* globe = hittable_sphere_new(world, (point3){0, 0, 0}, 2, (struct material*)earthsurface);
    return world;
}


///Perlin ball with simple light
hittable_list* simple_light(){
    hittable_list* world = hittable_list_new(1024);
    texture* pertext = texture_noise_init_scaled(4);
    material* mat1 = material_lambertian_new(pertext);
    hittable* s1 = hittable_sphere_new(world, (point3){0, -1000, 0}, 1000, (struct material*)mat1);
    hittable* s2 = hittable_sphere_new(world, (point3){0,  2, 0}, 2, (struct material*)mat1);
    material* difflight = material_light_new_from_color((color){4,4,4});
    hittable* l1 = hittable_rect_new(world, 3, 5, 1, 3, 0,0, -2,  XY, (struct material*)difflight);
    return world;
}


///Cornel box
hittable_list* cornell_box(){
    hittable_list* world = hittable_list_new(1024);
    material* red = material_lambertian_new_from_color((color){0.65, 0.05, 0.05});
    material* whi = material_lambertian_new_from_color((color){0.73, 0.73, 0.73});
    material* gre = material_lambertian_new_from_color((color){0.12, 0.45, 0.15});
    material* lig = material_light_new_from_color((color){15,15,15});
    hittable* b1 = hittable_rect_new(world, 0,0, 0, 555, 0, 555, 555, YZ, (struct material*)gre);
    hittable* b2 = hittable_rect_new(world, 0,0, 0, 555, 0, 555, 0, YZ, (struct material*)red);

    hittable* b3 = hittable_rect_new(world, 213, 343, 0,0, 227, 332, 554, XZ, (struct material*)lig);
    hittable* b4 = hittable_rect_new(world, 0, 555, 0,0, 0, 555, 0, XZ, (struct material*)whi);
    hittable* b5 = hittable_rect_new(world, 0, 555, 0,0, 0, 555, 555, XZ, (struct material*)whi);

    hittable* b6 = hittable_rect_new(world, 0, 555, 0, 555, 0,0, 555, XY, (struct material*)whi);
    hittable* box1 = hittable_box_new(NULL, (point3){0,0,0}, (point3){165,330,165}, (struct material*)whi);
    box1 = hittable_rotate_init(NULL, box1, 15, X);
    box1 = hittable_translate_init(world, box1, (vec3){265, 0, 295});
    hittable* box2 = hittable_box_new(NULL, (point3){0,0,0}, (point3){165,165,165}, (struct material*)whi);
    box2 = hittable_rotate_init(NULL, box2, -18, Y);
    box2 = hittable_translate_init(world, box2, (vec3){130, 0, 65});
    return world;
}


///Cornell box with smoke boxes
hittable_list* cornell_smoke(){
    hittable_list* world = hittable_list_new(1024);
    material* red = material_lambertian_new_from_color((color){0.65, 0.05, 0.05});
    material* whi = material_lambertian_new_from_color((color){0.73, 0.73, 0.73});
    material* gre = material_lambertian_new_from_color((color){0.12, 0.45, 0.15});
    material* lig = material_light_new_from_color((color){15,15,15});
    hittable* b1 = hittable_rect_new(world, 0,0, 0, 555, 0, 555, 555, YZ, (struct material*)gre);
    hittable* b2 = hittable_rect_new(world, 0,0, 0, 555, 0, 555, 0, YZ, (struct material*)red);
    hittable* b3 = hittable_rect_new(world, 213, 343, 0,0, 227, 332, 554, XZ, (struct material*)lig);
    hittable* b4 = hittable_rect_new(world, 0, 555, 0,0, 0, 555, 0, XZ, (struct material*)whi);
    hittable* b5 = hittable_rect_new(world, 0, 555, 0,0, 0, 555, 555, XZ, (struct material*)whi);
    hittable* b6 = hittable_rect_new(world, 0, 555, 0, 555, 0,0, 555, XY, (struct material*)whi);
    hittable* box1 = hittable_box_new(NULL, (point3){0,0,0}, (point3){165,330,165}, (struct material*)whi);
    box1 = hittable_rotate_init(NULL, box1, 15, Y);
    box1 = hittable_translate_init(NULL, box1, (vec3){265, 0, 295});
    box1 = hittable_constant_medium_init_c(world, box1, 0.01, (color){0,0,0});
    hittable* box2 = hittable_box_new(NULL, (point3){0,0,0}, (point3){165,165,165}, (struct material*)whi);
    box2 = hittable_rotate_init(NULL, box2, -18, Y);
    box2 = hittable_translate_init(NULL, box2, (vec3){130, 0, 65});
    box2 = hittable_constant_medium_init_c(world, box2, 0.01, (color){1,1,1});
    return world;
}


///Final scene with almost everything
hittable_list* final_scene(){
    //Setup ground with boxes
    hittable_list* boxes1 = hittable_list_new_no_gc(1024);
    material* ground = material_lambertian_new_from_color((color){0.48,0.83,0.53});
    const int boxes_per_side = 20;
    for (int i = 0; i < boxes_per_side; i++) {
        for (int j = 0; j < boxes_per_side; j++) {
            double w = 100.0;
            double x0 = -1000.0 + i*w;
            double z0 = -1000.0 + j*w;
            double y0 = 0.0;
            double x1 = x0 + w;
            double y1 = random_double_scaled(1,101);
            double z1 = z0 + w;
            hittable_box_new(boxes1, (point3){x0,y0,z0}, (point3){x1,y1,z1}, (struct material*)ground);
        }
    }

    //Wrap boxes into bvh tree
    hittable_list* objects = hittable_list_new(1024);
    hittable* node = bvh_node_init(objects, boxes1, 0, 1);

    //Setup light
    material* light = material_light_new_from_color((color){7,7,7});
    hittable_rect_new(objects, 123, 423, 0,0, 147, 412, 554, XZ, (struct material*)light);

    //Moving sphere
    point3 center1 = (point3){400, 400, 200};
    point3 center2 = vec3c_sum(center1, (vec3){30,0,0});
    material* moving_sphere_material = material_lambertian_new_from_color((color){0.7,0.3,0.1});
    hittable_moving_sphere_new(objects, center1, center2, 0, 1, 50, (struct material *)moving_sphere_material);

    //Glass ball
    hittable_sphere_new(objects, (point3){260,150,45}, 50, (struct material *)material_dielectric_new(1.5));

    //Metal ball
    hittable_sphere_new(objects, (point3){0,150,145}, 50, (struct material *)material_metal_new((color){0.8,0.8,0.8}, 1.0));

    //Gas ball
    hittable* boundary = hittable_sphere_new(NULL, (point3){360,150,145}, 70, (struct material *)material_dielectric_new(1.5));
    hittable_constant_medium_init_c(objects, boundary, 0.025, (color){0.2,0.4,0.9});

    //Huge light smoke area
    boundary = hittable_sphere_new(NULL, (point3){0,0,0}, 5000, (struct material *)material_dielectric_new(1.5));
    hittable_constant_medium_init_c(objects, boundary, 0.0001, (color){1,1,1});

    //Perline noise
    texture* pertext = texture_noise_init_scaled(5);
    hittable_sphere_new(objects, (point3){220,280,300}, 80, (struct material *)material_lambertian_new(pertext));

    //1000 small spheres grouped, translated, and rotated
    hittable_list* boxes2 = hittable_list_new_no_gc(1024);
    material* white = material_lambertian_new_from_color((color){0.73, 0.73, 0.73});
    int ns = 1000;
    for (int j = 0; j < ns; j++) {hittable_sphere_new(boxes2, vec3_random_scaled(0,165), 10, (struct material *)white);}
    hittable* bvh = bvh_node_init(NULL, boxes2, 0, 1);
    hittable* rotated = hittable_rotate_init(NULL, bvh, 15, Y);
    hittable* translated = hittable_translate_init(objects, rotated, (vec3){-100, 270, 395});
    return objects;
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

    //Define camera
    double ASPECT_RATIO = 16.0/9.0;
    double vfov = 40;
    point3 lookfrom;
    point3 lookat;
    vec3 vup = (vec3){0, 1, 0};
    double dist_to_focus = 10;
    double aperture = 0.0;
    color background = {0,0,0};

    //Define objects list
    hittable_list* world;
    switch(6){
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

        case 6:
            world = cornell_box();
            ASPECT_RATIO = 1.0;
            background = (color){0.0, 0.0, 0.0};
            lookfrom = (point3){278, 278, -800};
            lookat   = (point3){278, 278, 0};
            vfov = 40.0;
            break;

        case 7:
            world = cornell_smoke();
            ASPECT_RATIO = 1.0;
            lookfrom = (point3){278, 278, -800};
            lookat   = (point3){278, 278, 0};
            vfov = 40.0;
            break;

        default:
        case 8:
            world = final_scene();
            ASPECT_RATIO = 1.0;
            background = (color){0.0, 0.0, 0.0};
            lookfrom = (point3){478, 278, -600};
            lookat   = (point3){278, 278, 0};
            vfov = 40.0;
            break;
    }

    //Init camera
    const int IMAGE_WIDTH = 1024/8;
    const int IMAGE_HEIGHT = (int)(IMAGE_WIDTH / ASPECT_RATIO);
    const int samples_per_pixel = 32; //how many ray per pixel
    const int max_recursion_depth = 8; //how deep the ray scattering goes
    camera* c = camera_new(lookfrom, lookat, vup, vfov, ASPECT_RATIO, aperture, dist_to_focus, 0.0, 1.0);


    /**
     *
     * Rendering
     *
     **/

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
                    pixel_color = vec3c_sum(pixel_color, ray_color(&r, &background, world, max_recursion_depth));
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

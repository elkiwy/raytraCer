#define PROGRAM_FILE "src/program.cl"
#define ARRAY_SIZE 1024

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/time.h>

#include <OpenCL/cl.h>


#include "utils.h"
#include "camera.h"
#include "objects.h"
#include "renderer.h"










int main() {
    Material* mat_red = make_material_lambertian(make_texture_solid(0.65, 0.05, 0.05));
    Material* mat_white = make_material_lambertian(make_texture_solid(0.73, 0.73, 0.73));
    Material* mat_green = make_material_lambertian(make_texture_solid(0.12, 0.45, 0.15));
    Material* mat_light = make_material_light(make_texture_solid(15, 15, 15));
    //Material* mat_isotropic = make_material_isotropic(make_texture_solid(1, 1, 1));
    //Material* mat_glass = make_material_glass(1.5);
    //Material* mat_metal = make_material_metal(0, make_texture_solid(0.8, 0.85, 0.88));

    Object* r1 = make_rect(0,0, 555,555, 555, AXIS_YZ, mat_green);
    Object* r2 = make_rect(0,0, 555,555,   0, AXIS_YZ, mat_red);
    Object* light1 = make_rect(213,227, 343,332, 554, AXIS_XZ, mat_light);
    Object* light2 = make_flip_face(light1);
    Object* r3 = make_rect(0,0, 555,555,   0, AXIS_XZ, mat_white);
    Object* r4 = make_rect(0,0, 555,555, 555, AXIS_XZ, mat_white);
    Object* r5 = make_rect(0,0, 555,555, 555, AXIS_XY, mat_white);
    Object* box1 = make_box(0,0,0,  165,330,165, mat_white);
    box1 = make_rotated(box1, AXIS_Y, 15);
    box1 = make_translated(box1, 265,0,295);
    //cl_float16 box2 = make_box(wrapped_objs, &wrapped_obj_ind, 0,0,0,  165,165,165, mat_white);
    //box2 = make_rotated(wrapped_objs, &wrapped_obj_ind, box2, AXIS_Y, -18);
    //box2 = make_translated(wrapped_objs, &wrapped_obj_ind, box2, 130,0,65);
    //add_object(box2, objs, &obj_ind);
    //add_object(make_constant_medium_sphere(    275, 300, 275, 2500,  0.0005, mat_isotropic), objs, &obj_ind);

    ObjectList* objectList = make_objectList(16);
    add_object_to_objectList(r1, objectList);
    add_object_to_objectList(r2, objectList);
    add_object_to_objectList(light2, objectList);
    add_object_to_objectList(r3, objectList);
    add_object_to_objectList(r4, objectList);
    add_object_to_objectList(r5, objectList);
    add_object_to_objectList(box1, objectList);

    ObjectList* lightList  = make_objectList(16);
    add_object_to_objectList(light1, lightList);





    scene* s = make_scene(objectList, lightList, 128);

    scene_settings settings;
    settings.res_width         = 256;
    settings.res_height        = settings.res_width;
    settings.samples_per_pixel = 128;
    settings.iterations        = 8;

    float3 from = {278, 278, -800};
    float3   to = {278, 278, 0};
    float3  vup = { 0, 1, 0};
    float dist_to_focus = 10.0;
    float aperture = 0.0;
    float vfov = 40.0;
    camera* cam = make_camera(from, to, vup, vfov, 1.0, aperture, dist_to_focus);

    render_scene(s, cam, &settings);

    free_camera(cam);
    free_scene(s);

    return 0;
}

#ifndef __OBJECTS_H_
#define __OBJECTS_H_



#include <OpenCL/cl.h>




/*
Materials:
    0 -> type_flag
    1,2,3 -> albedo
    4 -> fuzz
    5 -> ir
    6 -> texture_index
*/

/*
Objects / Lights:
    0 -> type_flag
    1,2,3 -> point1
    4,5,6 -> point2
    7 -> radius / axis
    8 -> material_index
    9 -> rect_k / rotation_angle / density

    10 -> obj_ptr1
    11 -> obj_ptr2
    12 -> obj_ptr3
    13 -> obj_ptr4
    14 -> obj_ptr5
    15 -> obj_ptr6
*/

/*
Texture:
    0 -> type_flag
    1,2,3 -> color
    4 -> scale
*/


typedef enum{
    MAT_LAMBERTIAN = 0,
    MAT_METAL = 1,
    MAT_DIELECTRIC = 2,
    MAT_LIGHT = 3,
    MAT_ISOTROPIC = 4
} material_type;

typedef enum{
    OBJ_SPHERE = 0,
    OBJ_RECT = 1,
    OBJ_BOX = 2,
    OBJ_ROTATED = 3,
    OBJ_TRANSLATED = 4,
    OBJ_CONSTANT_MEDIUM = 5,
    OBJ_FLIP_FACE = 6
} object_type;

typedef enum{
    TEX_SOLID = 0,
    TEX_PERLIN = 1
} texture_type;


typedef enum{AXIS_XY = 0, AXIS_YZ = 1, AXIS_XZ = 2} axis;
typedef enum{AXIS_X = 0, AXIS_Y = 1, AXIS_Z = 2} rot_axis;





typedef struct Texture{
    cl_float16 data;
}Texture;


typedef struct Material{
    cl_float16 data;
    Texture* texture;
}Material;


typedef struct Object{
    cl_float16 data;
    Material* material;
    void* wrapped_objs[6];
}Object;





typedef struct ObjectList{
    Object** objects;
    int used;
    int size;
}ObjectList;


ObjectList* make_objectList(int size);
void add_object_to_objectList(Object* obj, ObjectList* list);











int add_to_list(cl_float16 item, cl_float16* arr, int* ind);
void free_object(Object* obj);



/**
 *
 * Textures creation
 *
 * */

Texture* make_texture_solid(float r,float g,float b);



/**
 *
 * Materials creation
 *
 * */

Material* make_material_lambertian(Texture* texture);
Material* make_material_light(Texture* texture);
Material* make_material_isotropic(Texture* texture);
Material* make_material_metal(float fuzz, Texture* texture);
Material* make_material_dielectric(float ir);



/**
 *
 * Object creation
 *
 * */

Object* make_sphere(float cx,float cy,float cz, float r, Material* mat);
Object* make_rect(float c00,float c01, float c10,float c11, float z, axis a, Material* mat);
Object* make_box(float c0x,float c0y, float c0z,float c1x, float c1y, float c1z, Material* mat);
Object* make_rotated(Object* obj, rot_axis ax, float angle);
Object* make_translated(Object* obj, float ox, float oy, float oz);
Object* make_flip_face(Object* obj);
Object* make_constant_medium_sphere(float cx, float cy, float cz, float r, float density, Material* mat);



/**
 * Objects utilities
 * */

int sameClFloat16(cl_float16* a, cl_float16* b);
int addItemIgnoringDuplicates(cl_float16 item, cl_float16* arr, int* index);
void addToUniques(Object* o, Object** out_objs, Material** out_mats, Texture** out_texs, int* out_objs_count, int* out_mats_count, int* out_texs_count);
int packObjectToGPUArrays(Object* obj, cl_float16* objs, int* objs_ind, cl_float16* mats, int* mats_ind, cl_float16* texs, int* texs_ind, cl_float16* wraps, int* wraps_ind);

#endif // __OBJECTS_H_

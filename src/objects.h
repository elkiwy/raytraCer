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




int add_to_list(cl_float16 item, cl_float16* arr, int* ind){
    arr[*ind] = item;
    *ind = *ind + 1;
    return *ind - 1;
}



/**
 *
 * Textures creation
 *
 * */

cl_float16 make_texture_solid(float r,float g,float b){
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = TEX_SOLID;
    o.s[1] = r;
    o.s[2] = g;
    o.s[3] = b;
    return o;
}






/**
 *
 * Materials creation
 *
 * */

cl_float16 make_material_lambertian(int texture_id){
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = MAT_LAMBERTIAN;
    o.s[6] = texture_id;
    return o;
}



cl_float16 make_material_light(int texture_id){
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = MAT_LIGHT;
    o.s[6] = texture_id;
    return o;
}



cl_float16 make_material_isotropic(int texture_id){
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = MAT_ISOTROPIC;
    o.s[6] = texture_id;
    return o;
}



cl_float16 make_material_metal(float fuzz, int texture_id){
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = MAT_ISOTROPIC;
    o.s[4] = fuzz;
    o.s[6] = texture_id;
    return o;
}



cl_float16 make_material_glass(float ir){
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = MAT_DIELECTRIC;
    o.s[5] = ir;
    return o;
}











/**
 *
 * Object creation
 *
 * */




cl_float16 make_sphere(float cx,float cy,float cz, float r, int mat_ind){
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = OBJ_SPHERE;
    o.s[1] = cx;
    o.s[2] = cy;
    o.s[3] = cz;
    o.s[7] = r;
    o.s[8] = mat_ind;
    return o;
}



cl_float16 make_rect(float c00,float c01, float c10,float c11, float z, axis a, int mat_ind){
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = OBJ_RECT;
    int c00_pos, c01_pos, c10_pos, c11_pos;
    if (a==AXIS_XY){c00_pos = 1; c01_pos = 2; c10_pos = 4; c11_pos = 5;
    }else if (a==AXIS_YZ){c00_pos = 2; c01_pos = 3; c10_pos = 5; c11_pos = 6;
    }else{c00_pos = 1; c01_pos = 3; c10_pos = 4; c11_pos = 6;}
    o.s[c00_pos] = c00;
    o.s[c01_pos] = c01;
    o.s[c10_pos] = c10;
    o.s[c11_pos] = c11;
    o.s[9] = z;
    o.s[7] = a;
    o.s[8] = mat_ind;
    return o;
}



cl_float16 make_box(cl_float16* wrapped_objs, int* wrapped_obj_ind, float c0x,float c0y, float c0z,float c1x, float c1y, float c1z, int mat_ind){
    int ptr1 = add_to_list(make_rect(c0x,c0y, c1x,c1y, c1z, AXIS_XY, mat_ind), wrapped_objs, wrapped_obj_ind);
    int ptr2 = add_to_list(make_rect(c0x,c0y, c1x,c1y, c0z, AXIS_XY, mat_ind), wrapped_objs, wrapped_obj_ind);
    int ptr3 = add_to_list(make_rect(c0x,c0z, c1x,c1z, c1y, AXIS_XZ, mat_ind), wrapped_objs, wrapped_obj_ind);
    int ptr4 = add_to_list(make_rect(c0x,c0z, c1x,c1z, c0y, AXIS_XZ, mat_ind), wrapped_objs, wrapped_obj_ind);
    int ptr5 = add_to_list(make_rect(c0y,c0z, c1y,c1z, c1x, AXIS_YZ, mat_ind), wrapped_objs, wrapped_obj_ind);
    int ptr6 = add_to_list(make_rect(c0y,c0z, c1y,c1z, c0x, AXIS_YZ, mat_ind), wrapped_objs, wrapped_obj_ind);

    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = OBJ_BOX;
    o.s[1] = c0x;
    o.s[2] = c0y;
    o.s[3] = c0z;
    o.s[4] = c1x;
    o.s[5] = c1y;
    o.s[6] = c1z;
    o.s[8] = mat_ind;
    o.s[10] = ptr1;
    o.s[11] = ptr2;
    o.s[12] = ptr3;
    o.s[13] = ptr4;
    o.s[14] = ptr5;
    o.s[15] = ptr6;
    return o;
}



cl_float16 make_rotated(cl_float16* wrapped_objs, int* wrapped_obj_ind, cl_float16 obj, rot_axis ax, float angle){
    int ptr = add_to_list(obj, wrapped_objs, wrapped_obj_ind);
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = OBJ_ROTATED;
    o.s[7] = ax;
    o.s[9] = angle;
    o.s[10] = ptr;
    return o;
}



cl_float16 make_translated(cl_float16* wrapped_objs, int* wrapped_obj_ind, cl_float16 obj, float ox, float oy, float oz){
    int ptr = add_to_list(obj, wrapped_objs, wrapped_obj_ind);
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = OBJ_TRANSLATED;
    o.s[1] = ox;
    o.s[2] = oy;
    o.s[3] = oz;
    o.s[10] = ptr;
    return o;
}



cl_float16 make_flip_face(cl_float16* wrapped_objs, int* wrapped_obj_ind, cl_float16 obj){
    int ptr = add_to_list(obj, wrapped_objs, wrapped_obj_ind);
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = OBJ_FLIP_FACE;
    o.s[10] = ptr;
    return o;
}



cl_float16 make_constant_medium_sphere(float cx, float cy, float cz, float r, float density, int mat_ind){
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = OBJ_CONSTANT_MEDIUM;
    o.s[1] = cx;
    o.s[2] = cy;
    o.s[3] = cz;
    o.s[7] = r;
    o.s[8] = mat_ind;
    o.s[9] = -1.0f/density;
    return o;
}


#endif // __OBJECTS_H_

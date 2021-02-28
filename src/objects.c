#include "objects.h"


ObjectList* make_objectList(int size){
    ObjectList* l = malloc(sizeof(ObjectList));
    l->objects = malloc(sizeof(Object*)*size);
    l->used = 0;
    l->size = size;
    return l;
}

void add_object_to_objectList(Object* obj, ObjectList* list){
    list->objects[list->used] = obj;
    list->used++;
}




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

Texture* make_texture_solid(float r,float g,float b){
    Texture* t = malloc(sizeof(Texture));
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = TEX_SOLID;
    o.s[1] = r;
    o.s[2] = g;
    o.s[3] = b;
    t->data = o;
    return t;
}






/**
 *
 * Materials creation
 *
 * */

Material* make_material_lambertian(Texture* texture){
    Material* m = malloc(sizeof(Material));
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = MAT_LAMBERTIAN;
    //o.s[6] = texture_id;
    m->data = o;
    m->texture = texture;
    return m;
}



Material* make_material_light(Texture* texture){
    Material* m = malloc(sizeof(Material));
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = MAT_LIGHT;
    //o.s[6] = texture_id;
    m->data = o;
    m->texture = texture;
    return m;
}



Material* make_material_isotropic(Texture* texture){
    Material* m = malloc(sizeof(Material));
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = MAT_ISOTROPIC;
    //o.s[6] = texture_id;
    m->data = o;
    m->texture = texture;
    return m;
}



Material* make_material_metal(float fuzz, Texture* texture){
    Material* m = malloc(sizeof(Material));
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = MAT_ISOTROPIC;
    o.s[4] = fuzz;
    //o.s[6] = texture_id;
    m->data = o;
    m->texture = texture;
    return m;
}



Material* make_material_dielectric(float ir){
    Material* m = malloc(sizeof(Material));
    m->texture = 0;
    cl_float16 o = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    o.s[0] = MAT_DIELECTRIC;
    o.s[5] = ir;
    m->data = o;
    return m;
}











/**
 *
 * Object creation
 *
 * */


Object* new_object(){
    Object* o = malloc(sizeof(Object));
    o->wrapped_objs[0] = 0;
    o->wrapped_objs[1] = 0;
    o->wrapped_objs[2] = 0;
    o->wrapped_objs[3] = 0;
    o->wrapped_objs[4] = 0;
    o->wrapped_objs[5] = 0;

    o->material = 0;

    o->data = (cl_float16){{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    return o;
}



Object* make_sphere(float cx,float cy,float cz, float r, Material* mat){
    Object* obj = new_object();
    obj->data.s[0] = OBJ_SPHERE;
    obj->data.s[1] = cx;
    obj->data.s[2] = cy;
    obj->data.s[3] = cz;
    obj->data.s[7] = r;
    obj->material = mat;
    return obj;
}



Object* make_rect(float c00,float c01, float c10,float c11, float z, axis a, Material* mat){
    Object* obj = new_object();
    obj->data.s[0] = OBJ_RECT;
    int c00_pos, c01_pos, c10_pos, c11_pos;
    if (a==AXIS_XY){      c00_pos = 1; c01_pos = 2; c10_pos = 4; c11_pos = 5;
    }else if (a==AXIS_YZ){c00_pos = 2; c01_pos = 3; c10_pos = 5; c11_pos = 6;
    }else{                c00_pos = 1; c01_pos = 3; c10_pos = 4; c11_pos = 6;}
    obj->data.s[c00_pos] = c00;
    obj->data.s[c01_pos] = c01;
    obj->data.s[c10_pos] = c10;
    obj->data.s[c11_pos] = c11;
    obj->data.s[9] = z;
    obj->data.s[7] = a;
    obj->material = mat;
    return obj;
}



Object* make_box(float c0x,float c0y, float c0z,float c1x, float c1y, float c1z, Material* mat){
    Object* obj = new_object();
    Object* ptr1 = make_rect(c0x,c0y, c1x,c1y, c1z, AXIS_XY, mat);
    Object* ptr2 = make_rect(c0x,c0y, c1x,c1y, c0z, AXIS_XY, mat);
    Object* ptr3 = make_rect(c0x,c0z, c1x,c1z, c1y, AXIS_XZ, mat);
    Object* ptr4 = make_rect(c0x,c0z, c1x,c1z, c0y, AXIS_XZ, mat);
    Object* ptr5 = make_rect(c0y,c0z, c1y,c1z, c1x, AXIS_YZ, mat);
    Object* ptr6 = make_rect(c0y,c0z, c1y,c1z, c0x, AXIS_YZ, mat);

    obj->wrapped_objs[0] = ptr1;
    obj->wrapped_objs[1] = ptr2;
    obj->wrapped_objs[2] = ptr3;
    obj->wrapped_objs[3] = ptr4;
    obj->wrapped_objs[4] = ptr5;
    obj->wrapped_objs[5] = ptr6;

    obj->data.s[0] = OBJ_BOX;
    obj->data.s[1] = c0x;
    obj->data.s[2] = c0y;
    obj->data.s[3] = c0z;
    obj->data.s[4] = c1x;
    obj->data.s[5] = c1y;
    obj->data.s[6] = c1z;
    obj->material = mat;
    return obj;
}



Object* make_rotated(Object* obj_to_wrap, rot_axis ax, float angle){
    Object* obj = new_object();
    obj->wrapped_objs[0] = obj_to_wrap;
    obj->data.s[0] = OBJ_ROTATED;
    obj->data.s[7] = ax;
    obj->data.s[9] = angle;
    return obj;
}



Object* make_translated(Object* obj_to_wrap, float ox, float oy, float oz){
    Object* obj = new_object();
    obj->wrapped_objs[0] = obj_to_wrap;
    obj->data.s[0] = OBJ_TRANSLATED;
    obj->data.s[1] = ox;
    obj->data.s[2] = oy;
    obj->data.s[3] = oz;
    return obj;
}



Object* make_flip_face(Object* obj_to_wrap){
    Object* obj = new_object();
    obj->wrapped_objs[0] = obj_to_wrap;
    obj->data.s[0] = OBJ_FLIP_FACE;
    return obj;
}



Object* make_constant_medium_sphere(float cx, float cy, float cz, float r, float density, Material* mat){
    Object* obj = new_object();
    obj->data.s[0] = OBJ_CONSTANT_MEDIUM;
    obj->data.s[1] = cx;
    obj->data.s[2] = cy;
    obj->data.s[3] = cz;
    obj->data.s[7] = r;
    obj->data.s[9] = -1.0f/density;
    obj->material = mat;
    return obj;
}





/**
 * Objects utilities
 * */


int sameClFloat16(cl_float16* a, cl_float16* b){
    for(int i=0;i<16;i++){if (a->s[i] != b->s[i]){return 0;}}
    return 1;
}



int addItemIgnoringDuplicates(cl_float16 item, cl_float16* arr, int* index){
    //Search if is already in the array
    for(int i=0;i<*index;i++){if (sameClFloat16(&item, &arr[i])){return i;}}

    //Not found
    arr[*index] = item;
    *index = *index + 1;
    return *index - 1;
}



int packObjectToGPUArrays(Object* obj, cl_float16* objs, int* objs_ind, cl_float16* mats, int* mats_ind, cl_float16* texs, int* texs_ind, cl_float16* wraps, int* wraps_ind){
    if(obj->material != 0){
        Material* mat = obj->material;
        //Add texture
        if (mat->texture != 0){mat->data.s[6] = addItemIgnoringDuplicates(mat->texture->data, texs, texs_ind);}
        //Add Material
        obj->data.s[8] = addItemIgnoringDuplicates(mat->data, mats, mats_ind);
    }

    //Add wrapped objects if any
    for(int j=0; j<6; j++){
        if (obj->wrapped_objs[j] != 0){
            obj->data.s[10+j] = packObjectToGPUArrays(obj->wrapped_objs[j], NULL, NULL, mats, mats_ind, texs, texs_ind, wraps, wraps_ind);
        }
    }

    //Add Object
    int obj_id;
    if(objs == NULL){obj_id = addItemIgnoringDuplicates(obj->data, wraps, wraps_ind);
    }else{obj_id = addItemIgnoringDuplicates(obj->data, objs, objs_ind);}
    return obj_id;
}



void addToUniques(Object* o, Object** out_objs, Material** out_mats, Texture** out_texs, int* out_objs_count, int* out_mats_count, int* out_texs_count){
    if(o->material != 0){
        Material* m = o->material;

        //Texture
        if (m->texture != 0){
            Texture* t = m->texture;
            int already_found = 0;
            for(int j=0; j<*out_texs_count; j++){
                if (t == out_texs[j]){already_found = 1;}
            }
            if (!already_found){
                out_texs[*out_texs_count] = t;
                *out_texs_count = *out_texs_count + 1;
            }
        }

        //Material
        int already_found = 0;
        for(int j=0; j<*out_mats_count; j++){
            if (m == out_mats[j]){already_found = 1;}
        }
        if (!already_found){
            out_mats[*out_mats_count] = m;
            *out_mats_count = *out_mats_count + 1;
        }
    }

    //Object
    int already_found = 0;
    for(int j=0; j<*out_objs_count; j++){
        if (o == out_objs[j]){already_found = 1;}
    }
    if (!already_found){
        out_objs[*out_objs_count] = o;
        *out_objs_count = *out_objs_count + 1;

        for(int j=0; j<6; j++){
            if (o->wrapped_objs[j] != 0){
                addToUniques(o->wrapped_objs[j], out_objs, out_mats, out_texs, out_objs_count, out_mats_count, out_texs_count);
            }
        }
    }
}

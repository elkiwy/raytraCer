#include "hittable.h"

#include "hittable_list.h"


/**
 *
 * Hit record
 *
 * */

///Set the front_face and normal property of the hit record object based on the ray and normal
inline static void hit_record_set_facenormal(hit_record* rec, ray* r, vec3* outward_normal){
    rec->front_face = vec3_dot(&(r->dir), outward_normal) < 0;
    rec->normal = rec->front_face ? vec3_copy(outward_normal) : vec3_flip(outward_normal);
}

///Print the hit record object
void hit_record_print(hit_record* h){
    printf("Rec: normal: [%.2f, %.2f, %.2f] p: [%.2f, %.2f, %.2f] t: %.2f, front_face: %d\n", h->normal.x, h->normal.y, h->normal.z, h->p.x, h->p.y, h->p.z, h->t, h->front_face );
}








/**
 *
 * Utility
 *
 * */

///Rotate a point in space, useful on hittable_rotate
void rotatePoint(vec3 original, double sin_theta, double cos_theta, rotation_axis axis, vec3* output){
    if (axis == X){
        output->x =  cos_theta * original.x  -  sin_theta * original.y;
        output->y =  sin_theta * original.x  +  cos_theta * original.y;
    }else if (axis == Y){
        output->x =  cos_theta * original.x  +  sin_theta * original.z;
        output->z = -sin_theta * original.x  +  cos_theta * original.z;
    }else{
        output->y =  cos_theta * original.y  -  sin_theta * original.z;
        output->z =  sin_theta * original.y  +  cos_theta * original.z;
    }
}



/**
 *
 * BVH AABB sorting
 *
 * */

int hittable_bounding_box(hittable* a, double t0, double t1, aabb* output_box);
inline static int box_compare(hittable* a, hittable* b, int axis){
    aabb box_a, box_b;
    if(!hittable_bounding_box(a, 0, 0, &box_a) || !hittable_bounding_box(b, 0, 0, &box_b)){
        printf("Failed to get bounding box in bvh_node comparator.\n"); fflush(stdout);
        return 0;}
    if (axis==0){       return box_a.minimum.x < box_b.minimum.x;
    }else if (axis==1){ return box_a.minimum.y < box_b.minimum.y;
    }else{              return box_a.minimum.z < box_b.minimum.z;}}
int box_x_compare(const void* a, const void* b){return box_compare(*((hittable**)a), *((hittable**)b), 0);}
int box_y_compare(const void* a, const void* b){return box_compare(*((hittable**)a), *((hittable**)b), 1);}
int box_z_compare(const void* a, const void* b){return box_compare(*((hittable**)a), *((hittable**)b), 2);}





/**
 *
 * Hittable Objects allocators
 *
 * */


/**
 * Objects
 * */

///Create a generic hittable wrapper object
hittable* hittable_generic_new(hittable_list* world, hittable_type t, void* obj){
    hittable* h = malloc(sizeof(hittable));
    h->t = t; h->obj = obj;
    hittable_list_add(world, h);
    return h;
}


///Create a new sphere object
hittable* hittable_sphere_new(struct hittable_list* world, point3 center, double r, struct material* mat){
    //Create the sphere
    sphere* s = malloc(sizeof(sphere));
    s->center = center;
    s->r = r;
    s->mat = mat;

    //Make generic wrapper and return
    return hittable_generic_new(world, HITTABLE_SPHERE, s);
}


///Create a new moving sphere object
hittable* hittable_moving_sphere_new(struct hittable_list* world, point3 center0, point3 center1, double t0, double t1, double r, struct material* mat){
    //Create the sphere
    moving_sphere* s = malloc(sizeof(moving_sphere));
    s->center0 = center0;
    s->center1 = center1;
    s->time0 = t0;
    s->time1 = t1;
    s->r = r;
    s->mat = mat;

    //Make generic wrapper and return
    return hittable_generic_new(world, HITTABLE_MOVING_SPHERE, s);
}


///Create a XY axis aligned rect
hittable* hittable_rect_new(struct hittable_list* world, double x0, double x1, double y0, double y1, double z0, double z1, double k, rect_axis axis, struct material* mat){
    rect* r = malloc(sizeof(rect));
    r->x0 = x0; r->x1 = x1;
    r->y0 = y0; r->y1 = y1;
    r->z0 = z0; r->z1 = z1;
    r->k = k;
    r->axis = axis;
    r->mat = mat;
    //Make generic wrapper and return
    return hittable_generic_new(world, HITTABLE_RECT, r);
}


///Create a 3d box object
hittable* hittable_box_new(struct hittable_list* world, point3 min, point3 max, struct material* mat){
    box* b = malloc(sizeof(box));
    b->box_max = max;
    b->box_min = min;

    //Make all 6 sides
    b->sides = (struct hittable_list*)hittable_list_new(6);
    hittable_rect_new(b->sides, min.x, max.x, min.y, max.y, 0, 0, max.z, XY, mat);
    hittable_rect_new(b->sides, min.x, max.x, min.y, max.y, 0, 0, min.z, XY, mat);
    hittable_rect_new(b->sides, min.x, max.x, 0, 0, min.z, max.z, max.y, XZ, mat);
    hittable_rect_new(b->sides, min.x, max.x, 0, 0, min.z, max.z, min.y, XZ, mat);
    hittable_rect_new(b->sides, 0, 0, min.y, max.y, min.z, max.z, max.x, YZ, mat);
    hittable_rect_new(b->sides, 0, 0, min.y, max.y, min.z, max.z, min.x, YZ, mat);

    //Make generic wrapper and return
    return hittable_generic_new(world, HITTABLE_BOX, b);
}


/**
 * BVH Nodes
 * */

///Create a BVH node recursevly selecting src objects
hittable* bvh_node_constructor(hittable** src_objects, int src_objects_count, int start, int end, double t0, double t1){
    //Make node instance
    bvh_node* node = malloc(sizeof(bvh_node));

    //Create copy of src_object to be able to sort it without messing with the original order
    hittable** objects = malloc(sizeof(hittable*)*src_objects_count);
    memcpy(objects, src_objects, sizeof(hittable*)*src_objects_count);

    //Select a random axis
    int raxis = random_double();
    int(*comparator)(const void*, const void*);
    if (raxis < 0.33){comparator = &box_x_compare;
    }else if (raxis < 0.66){comparator = &box_y_compare;
    }else{comparator = &box_z_compare;}

    //Create leafs based on how many items i have to organize
    int object_span = end-start;
    if(object_span == 1){
        //Only one: make both leafs the same object
        node->left = objects[start];
        node->right = node->left;
    }else if (object_span == 2){
        //Two items, sort those and assign
        int ind1, ind2;
        if(comparator(&objects[start], &objects[start+1])){ind1 = start; ind2 = start+1;
        }else{ind1 = start+1; ind2 = start;}
        node->left = objects[ind1];
        node->right = objects[ind2];
    }else{
        //Multiple items: recursively sort, divide et impera
        qsort(objects, src_objects_count, sizeof(objects[0]), comparator);
        int mid = start + object_span/2;
        node->left  = bvh_node_constructor(src_objects, src_objects_count, start, mid, t0, t1);
        node->right = bvh_node_constructor(src_objects, src_objects_count, mid, end, t0, t1);
    }

    //Make the two leafs AABB
    aabb box_left, box_right;
    if(!hittable_bounding_box(node->left, t0, t1, &box_left) || !hittable_bounding_box(node->right, t0, t1, &box_right)){
        printf("Failed to get bounding box in bvh_node constructor.\n"); fflush(stdout); return NULL;}

    //Make the current AABB
    node->box = surrounding_box(box_left, box_right);

    //Free the temporary object copy
    free(objects);

    //Make generic and return (NOTE: we do not add it to the current hittable_list!)
    hittable* h = malloc(sizeof(hittable));
    h->t = HITTABLE_BVH_NODE;
    h->obj = node;
    return h;
}


///Create the root of a BVH tree
hittable* bvh_node_init(struct hittable_list* world, struct hittable_list* objects, double t0, double t1){
    void* objs = hittable_list_objects(objects);
    int index = hittable_list_index(objects);
    hittable* o = bvh_node_constructor(objs, index, 0, index, t0, t1);
    hittable_list_add(world, o);
    return o;
}





/**
 * Transformation wrappers
 * */

///Create a translate wrapper around a hittable object to move it by offset
hittable* hittable_translate_init(struct hittable_list* world, hittable* obj, vec3 offset){
    translate* translated = malloc(sizeof(translate));
    translated->obj = obj;
    translated->offset = offset;
    //Make generic wrapper and return
    return hittable_generic_new(world, HITTABLE_TRANSLATE, translated);
}


///Create a rotate wrapper around a hittable object to rotate it
hittable* hittable_rotate_init(struct hittable_list* world, hittable* obj, double angle, rotation_axis axis){
    rotate* rotated = malloc(sizeof(rotate));
    rotated->obj = obj;
    rotated->axis = axis;

    //Save some handy parameters
    double radians = deg2rad(angle);
    rotated->sin_theta = sin(radians);
    rotated->cos_theta = cos(radians);
    rotated->hasbox = hittable_bounding_box(obj, 0, 1, &rotated->bbox);

    //Create the new AABB
    point3 min = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
    point3 max = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            for(int k=0; k<2; k++){
                double x = i * rotated->bbox.maximum.x + (1-i) * rotated->bbox.minimum.x;
                double y = j * rotated->bbox.maximum.y + (1-j) * rotated->bbox.minimum.y;
                double z = k * rotated->bbox.maximum.z + (1-k) * rotated->bbox.minimum.z;

                vec3 p = {x, y, z};
                vec3 tester = p;
                rotatePoint(p, rotated->sin_theta, rotated->cos_theta, rotated->axis, &tester);

                //double newx = rotated->cos_theta * x + rotated->sin_theta * z;
                //double newz = -rotated->sin_theta * x + rotated->cos_theta * z;

                //vec3 tester = {newx, y, newz};
                min.x = fmin(min.x, tester.x);
                max.x = fmax(max.x, tester.x);
                min.y = fmin(min.y, tester.y);
                max.y = fmax(max.y, tester.y);
                min.z = fmin(min.z, tester.z);
                max.z = fmax(max.z, tester.z);
            }
        }
    }
    rotated->bbox = (aabb){min, max};

    //Make generic wrapper and return
    return hittable_generic_new(world, HITTABLE_ROTATE, rotated);
}


///Create a constant medium wrapper for effect like smoke
hittable* hittable_constant_medium_init(struct hittable_list* world, hittable* b, double d, texture* a){
    constant_medium* m = malloc(sizeof(constant_medium));
    m->boundary = b;
    m->neg_inv_density = -1/d;
    m->phase_function = (struct material*)material_isotropic_new(a);
    //Make generic wrapper and return
    return hittable_generic_new(world, HITTABLE_CONSTANT_MEDIUM, m);}
///Create a constant medium wrapper for effect like smoke from color
hittable* hittable_constant_medium_init_c(struct hittable_list* world, hittable* b, double d, color c){
    texture* t = texture_solid_color_init(c);
    return hittable_constant_medium_init(world, b, d, t);
}


///Create a flip face wrapper around a hittable object
hittable* hittable_flip_face_init(struct hittable_list* world, hittable* obj){
    flip_face* flipped = malloc(sizeof(flip_face));
    flipped->obj = obj;
    //Make generic wrapper and return
    return hittable_generic_new(world, HITTABLE_FLIP_FACE, flipped);
}
















/**
 *
 * Hittable Objects deallocators
 *
 * */

///Free hittable object
void hittable_free(hittable* o){
    if(o->t == HITTABLE_BVH_NODE){
        bvh_node* n = o->obj;
        if (n->left == n->right){
            //Same leaf, free once
            hittable_free(n->left);
        }else{
            //Free both leafs
            hittable_free(n->left);
            hittable_free(n->right);
        }

    }else if(o->t == HITTABLE_TRANSLATE){
        translate* r = o->obj;
        hittable_free(r->obj);

    }else if(o->t == HITTABLE_ROTATE){
        rotate* r = o->obj;
        hittable_free(r->obj);

    }else if(o->t == HITTABLE_CONSTANT_MEDIUM){
        constant_medium* m = o->obj;
        hittable_free(m->boundary);

    }else if(o->t == HITTABLE_FLIP_FACE){
        flip_face* m = o->obj;
        hittable_free(m->obj);
    }
    free(o->obj);
    free(o);
}








/**
 *
 * UV mapping
 *
 * */

void sphere_get_uv(point3* p, double* u, double* v){
    double theta = acos(-p->y);
    double phi = atan2(-p->z, p->x) + PI;
    *u = phi / (2*PI);
    *v = theta / PI;
}








/**
 *
 * Hittable Objects hit testing
 *
 * */

//Forward declaration for cyclic dependencies
int hittable_hit(hittable* o, ray* r, double tmin, double tmax, hit_record* rec);


///Test intersection between line and spehere
int line_sphere_intersection(vec3 oc, vec3 dir, double rad, double tmin, double tmax, double* root){
    double a = vec3_length_squared(&dir);
    double half_b = vec3_dot(&oc, &dir);
    double c = vec3_length_squared(&oc) - rad * rad;
    double discriminant = half_b*half_b - a*c;

    //If discriminant is < 0, then there is no intersection...
    if (discriminant < 0){return 0;}
    //...else the sphere and line intersect, so we have to check if the intersection are between tmin and tmax

    //Find nearest root that lies in the acceptable range.
    double sqrtd = sqrt(discriminant);
    double tmp_root = (-half_b - sqrtd) / a;
    if (tmp_root < tmin || tmp_root > tmax){
        //Outside bounds, check the other root
        tmp_root = (-half_b + sqrtd) / a;
        if (tmp_root < tmin || tmp_root > tmax){
            return 0;
        }
    }

    *root = tmp_root;
    return 1;
}


///Calculate sphere/ray hit and return a hit_record and a success state (0/1)
int sphere_hit(sphere* s, ray* r, double tmin, double tmax, hit_record* rec){
    //Find the discriminant from the formula of intersection between a 3d line and a sphere
    vec3 oc = vec3_sub(&(r->orig), &(s->center));

    double root;
    if (line_sphere_intersection(oc, r->dir, s->r, tmin, tmax, &root)) {
        //Here we should have valid root, so we fill the data of the hit record
        rec->t = root;
        rec->p = ray_at(r, rec->t);
        vec3 outward_normal = vec3c_div_k(vec3_sub(&rec->p, &s->center), s->r);
        hit_record_set_facenormal(rec, r, &outward_normal);
        rec->mat = s->mat;
        sphere_get_uv(&outward_normal, &rec->u, &rec->v);
        return 1;
    }else{
        //Intersection not found
        return 0;
    }
}


///Fine the center of the sphere in a certain time t
point3 moving_sphere_center(moving_sphere* s, double t){
    double k = ((t - s->time0) / (s->time1 - s->time0));
    return vec3c_sum(vec3c_mul_k(vec3_sub(&s->center1, &s->center0), k), s->center0);
}


///Calculate sphere/ray hit and return a hit_record and a success state (0/1)
int moving_sphere_hit(moving_sphere* s, ray* r, double tmin, double tmax, hit_record* rec){
    //Find the discriminant from the formula of intersection between a 3d line and a sphere
    point3 center_in_time = moving_sphere_center(s, r->time);
    vec3 oc = vec3c_sub(r->orig, center_in_time);

    double root;
    if (line_sphere_intersection(oc, r->dir, s->r, tmin, tmax, &root)) {
        //Here we should have valid root, so we fill the data of the hit record
        rec->t = root;
        rec->p = ray_at(r, rec->t);
        vec3 outward_normal = vec3c_div_k(vec3_sub(&rec->p, &center_in_time), s->r);
        hit_record_set_facenormal(rec, r, &outward_normal);
        rec->mat = s->mat;
        sphere_get_uv(&outward_normal, &rec->u, &rec->v);
        return 1;
    }else{
        //Intersection not found
        return 0;
    }
}


///Test hit on a BVH node
int bvh_node_hit(bvh_node* node, ray* r, double tmin, double tmax, hit_record* rec){
    if(aabb_hit(&node->box, r, tmin, tmax) == 0){return 0;}
    int hit_left = hittable_hit(node->left, r, tmin, tmax, rec);
    int hit_right = hittable_hit(node->right, r, tmin, hit_left ? rec->t : tmax, rec);
    return hit_left || hit_right;
}


///Test hit on a axis aligned rect
int rect_hit(rect* rect, ray* r, double tmin, double tmax, hit_record* rec){
    //Check for plane ray intersection
    double t;
    if(rect->axis == XY){t = (rect->k - r->orig.z) / r->dir.z;
    }else if(rect->axis == XZ){t = (rect->k - r->orig.y) / r->dir.y;
    }else{t = (rect->k - r->orig.x) / r->dir.x;}
    if (t < tmin || t > tmax){return 0;}

    //Find if the point of intersection is in the boundary
    double x = r->orig.x + t * r->dir.x; double y = r->orig.y + t * r->dir.y; double z = r->orig.z + t * r->dir.z;
    if ((rect->axis == XZ && (x < rect->x0 || x> rect->x1 || z < rect->z0 || z > rect->z1))
     || (rect->axis == XY && (x < rect->x0 || x> rect->x1 || y < rect->y0 || y > rect->y1))
     || (rect->axis == YZ && (z < rect->z0 || z> rect->z1 || y < rect->y0 || y > rect->y1))) {return 0;}

    //Register the hit
    rec->u = rect->axis==YZ ? (y - rect->y0)/(rect->y1 - rect->y0) : (x - rect->x0)/(rect->x1 - rect->x0);
    rec->v = rect->axis==XY ? (y - rect->y0)/(rect->y1 - rect->y0) : (z - rect->z0)/(rect->z1 - rect->z0);
    rec->t = t;
    vec3 outward_normal = {rect->axis==YZ?1:0,rect->axis==XZ?1:0,rect->axis==XY?1:0};
    hit_record_set_facenormal(rec, r, &outward_normal);
    rec->mat = rect->mat;
    rec->p = ray_at(r, t);
    return 1;
}


///Test hit on a box
int box_hit(box* b, ray* r, double tmin, double tmax, hit_record* rec){
    return hittable_list_hit(b->sides, r, tmin, tmax, rec);
}


///Test hit on a translated object
int translate_hit(translate* b, ray* r, double tmin, double tmax, hit_record* rec){
    //Move the ray to simulate movement and test on the wrapped object
    ray moved_ray = {vec3c_sub(r->orig, b->offset), r->dir, r->time};
    if (!hittable_hit(b->obj, &moved_ray, tmin, tmax, rec)) return 0;

    //Move back the hit point
    rec->p = vec3c_sum(rec->p, b->offset);
    hit_record_set_facenormal(rec, &moved_ray, &rec->normal);
    return 1;
}


///Test hit on a rotated object
int rotate_hit(rotate* ob, ray* r, double tmin, double tmax, hit_record* rec){
    //Rotate the incoming ray
    point3 origin = r->orig;
    vec3 direction = r->dir;
    rotatePoint(r->orig, -ob->sin_theta, ob->cos_theta, ob->axis, &origin);
    rotatePoint(r->dir, -ob->sin_theta, ob->cos_theta, ob->axis, &direction);

    //Test on the wrapped object
    ray rotated_ray = {origin, direction, r->time};
    if (!hittable_hit(ob->obj, &rotated_ray, tmin, tmax, rec)) return 0;

    //Revert rotation on the normal
    vec3 p = rec->p;
    vec3 normal = rec->normal;
    rotatePoint(rec->p, ob->sin_theta, ob->cos_theta, ob->axis, &p);
    rotatePoint(rec->normal, ob->sin_theta, ob->cos_theta, ob->axis, &normal);
    rec->p = p;
    hit_record_set_facenormal(rec, &rotated_ray, &normal);
    return 1;
}


///Test hit on a constant medium
int constant_medium_hit(constant_medium* m, ray* r, double tmin, double tmax, hit_record* rec){
    //double Test hit
    hit_record rec1, rec2;
    if (!hittable_hit(m->boundary, r, -HUGE_VAL, HUGE_VAL, &rec1)){return 0;}
    if (!hittable_hit(m->boundary, r, rec1.t+0.0001, HUGE_VAL, &rec2)){return 0;}
    if (rec1.t < tmin) rec1.t = tmin;
    if (rec2.t > tmax) rec2.t = tmax;
    if (rec1.t >= rec2.t) return 0;
    if (rec1.t < 0){rec1.t = 0;}
    const double ray_length = vec3_length(&r->dir);
    const double distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
    const double hit_distance = m->neg_inv_density * log(random_double());
    if(hit_distance > distance_inside_boundary){return 0;}

    //Register hit
    rec->t = rec1.t + hit_distance / ray_length;
    rec->p = ray_at(r, rec->t);
    rec->normal = (vec3){1,0,0};
    rec->front_face = 1;
    rec->mat = m->phase_function;
    return 1;
}


int flip_face_hit(flip_face* f, ray* r, double tmin, double tmax, hit_record* rec){
    if (!hittable_hit(f->obj, r, tmin, tmax, rec))return 0;
    rec->front_face = !rec->front_face;
    return 1;
}


///Generic function to test for ray hit
int hittable_hit(hittable* o, ray* r, double tmin, double tmax, hit_record* rec){
    if (o->t == HITTABLE_SPHERE){
        return sphere_hit((sphere*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_MOVING_SPHERE){
        return moving_sphere_hit((moving_sphere*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_BVH_NODE){
        return bvh_node_hit((bvh_node*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_RECT){
        return rect_hit((rect*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_BOX){
        return box_hit((box*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_TRANSLATE){
        return translate_hit((translate*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_ROTATE){
        return rotate_hit((rotate*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_CONSTANT_MEDIUM){
        return constant_medium_hit((constant_medium*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_FLIP_FACE){
        return flip_face_hit((flip_face*)o->obj, r, tmin, tmax, rec);
    }else{printf("!! Not implemented hittable_hit for type %d\n", o->t);return 0;}
}














/**
 *
 * Hittable Objects Bounding box
 *
 * */


///Get sphere bounding box
int sphere_bounding_box(sphere* s, aabb* output_box){
    output_box->minimum = vec3c_sub(s->center, (vec3){s->r,s->r,s->r});
    output_box->maximum = vec3c_sum(s->center, (vec3){s->r,s->r,s->r});
    return 1;
}


///Get moving sphere bounding box
int moving_sphere_bounding_box(moving_sphere* s, double t0, double t1, aabb* output_box){
    vec3 area = {s->r, s->r, s->r};
    aabb box0 = {vec3c_sub(moving_sphere_center(s, t0), area), vec3c_sum(moving_sphere_center(s, t0), area)};
    aabb box1 = {vec3c_sub(moving_sphere_center(s, t1), area), vec3c_sum(moving_sphere_center(s, t1), area)};
    aabb res = surrounding_box(box0, box1);
    output_box->minimum = res.minimum;
    output_box->maximum = res.maximum;
    return 1;
}


///Get BVH bounding box
int bvh_node_bounding_box(bvh_node* node, aabb* output_box){
    *output_box = node->box;
    return 1;
}


///Get Rect bounding box
int rect_bounding_box(rect* rect, aabb* output_box){
    point3 a,b;
    if (rect->axis == XY){
        a = (point3){rect->x0, rect->y0, rect->k-0.0001};
        b = (point3){rect->x1, rect->y1, rect->k+0.0001};
    }else if(rect->axis == XZ){
        a = (point3){rect->x0, rect->k-0.0001, rect->z0};
        b = (point3){rect->x1, rect->k+0.0001, rect->z1};
    }else{
        a = (point3){rect->k-0.0001, rect->y0, rect->z0};
        b = (point3){rect->k+0.0001, rect->y1, rect->z1};
    }
    *output_box = (aabb){a, b};
    return 1;
}


///Get Box bounding box
int box_bounding_box(box* b, aabb* output_box){
    *output_box = (aabb){b->box_min, b->box_max};
    return 1;
}


///Get translated object bounding box
int translate_bounding_box(translate* t, double t0, double t1, aabb* output_box){
    if (!hittable_bounding_box(t->obj, t0, t1, output_box)) return 0;
    *output_box = (aabb){vec3c_sum(output_box->minimum, t->offset), vec3c_sum(output_box->maximum, t->offset)};
    return 1;
}


///Get rotated object bounding box
int rotate_bounding_box(rotate* ob, aabb* output_box){
    *output_box = ob->bbox;
    return ob->hasbox;
}


///Get const medium object bounding box
int constant_medium_bounding_box(constant_medium* m, double t0, double t1, aabb* output_box){
    return hittable_bounding_box(m->boundary, t0, t1, output_box);
}


///Get flip face object bounding box
int flip_face_bounding_box(flip_face* m, double t0, double t1, aabb* output_box){
    return hittable_bounding_box(m->obj, t0, t1, output_box);
}


///Generic function to test for ray hit
int hittable_bounding_box(hittable* o, double t0, double t1, aabb* output_box){
    if (o->t == HITTABLE_SPHERE){
        return sphere_bounding_box((sphere*)o->obj, output_box);
    }else if (o->t == HITTABLE_MOVING_SPHERE){
        return moving_sphere_bounding_box((moving_sphere*)o->obj, t0, t1, output_box);
    }else if (o->t == HITTABLE_BVH_NODE){
        return bvh_node_bounding_box((bvh_node*)o->obj, output_box);
    }else if (o->t == HITTABLE_RECT){
        return rect_bounding_box((rect*)o->obj, output_box);
    }else if (o->t == HITTABLE_BOX){
        return box_bounding_box((box*)o->obj, output_box);
    }else if (o->t == HITTABLE_TRANSLATE){
        return translate_bounding_box((translate*)o->obj, t0, t1, output_box);
    }else if (o->t == HITTABLE_ROTATE){
        return rotate_bounding_box((rotate*)o->obj, output_box);
    }else if (o->t == HITTABLE_CONSTANT_MEDIUM){
        return constant_medium_bounding_box((constant_medium*)o->obj, t0, t1, output_box);
    }else if (o->t == HITTABLE_FLIP_FACE){
        return flip_face_bounding_box((flip_face*)o->obj, t0, t1, output_box);
    }else{printf("!! Not implemented hittable_hit for type %d\n", o->t);return 0;}
}














/**
 *
 * Hittable Objects Probability density functions (PDFs)
 *
 * */


///PDF of a spehere
double hittable_sphere_pdf_value(sphere* s, point3 orig, vec3 v){
    hit_record rec;
    ray r = {orig, v, 0};
    if(!sphere_hit(s, &r, 0.001, HUGE_VAL, &rec)){return 0;}

    double cos_theta_max = sqrt(1 - s->r * s->r / vec3c_length_squared(vec3_sub(&s->center, &orig)));
    double solid_angle = 2 * PI * (1-cos_theta_max);
    return 1/solid_angle;
}


///PDF of a rect
double hittable_rect_pdf_value(rect* rect, point3 orig, vec3 v){
    hit_record rec;
    ray r = {orig, v, 0};
    if(!rect_hit(rect, &r, 0.0001, HUGE_VAL, &rec)){return 0;}
    double distance_squared = rec.t * rec.t * vec3_length_squared(&v);

    double area = 0;
    if(rect->axis == XZ){
        area = (rect->x1-rect->x0)*(rect->z1-rect->z0);
    }else if(rect->axis == XY){
        area = (rect->x1-rect->x0)*(rect->y1-rect->y0);
    }else if(rect->axis == YZ){
        area = (rect->y1-rect->y0)*(rect->z1-rect->z0);
    }
    double cosine = fabs(vec3_dot(&v, &rec.normal) / vec3_length(&v));
    return distance_squared / (cosine * area);
}


///Generic PDF value function call
double hittable_pdf_value(hittable* o, point3 orig, vec3 v){
    if (o->t == HITTABLE_SPHERE){
        return hittable_sphere_pdf_value((sphere*)o->obj, orig, v);
    }else if (o->t == HITTABLE_RECT){
        return hittable_rect_pdf_value((rect*)o->obj, orig, v);

    ///Wrappers
    }else if (o->t == HITTABLE_FLIP_FACE){
        return hittable_pdf_value(((flip_face*)o->obj)->obj, orig, v);

    }else{
        return 0.0;
    }
}








/**
 *
 * Hittable Objects Random point
 *
 * */

///Random point on a spehere
vec3 hittable_sphere_random(sphere* s, point3 orig){
    vec3 direction = vec3_sub(&s->center, &orig);
    double distance_squared = vec3_length_squared(&direction);
    onb uvw;
    onb_init_from_w(&uvw, &direction);

    //Random point on a sphere
    double r1 = random_double();
    double r2 = random_double();
    double z = 1 + r2 * (sqrt(1-s->r*s->r/distance_squared) - 1);
    double phi = 2 * PI * r1;
    double x = cos(phi)*sqrt(1-z*z);
    double y = sin(phi)*sqrt(1-z*z);
    return onb_localv(&uvw, (vec3){x, y, z});
}

///Random point on a rect
vec3 hittable_rect_random(rect* rect, point3 orig){
    point3 random_point;
    if(rect->axis == XZ){
        random_point = (point3){random_double_scaled(rect->x0, rect->x1), rect->k, random_double_scaled(rect->z0, rect->z1)};
    }else if(rect->axis == XY){
        random_point = (point3){random_double_scaled(rect->x0, rect->x1), random_double_scaled(rect->y0, rect->y1), rect->k};
    }else if(rect->axis == YZ){
        random_point = (point3){rect->k, random_double_scaled(rect->y0, rect->y1), random_double_scaled(rect->z0, rect->z1)};
    }
    return vec3c_sub(random_point, orig);
}

///Generic random point on an object
vec3 hittable_random(hittable* o, point3 orig){
    if (o->t == HITTABLE_SPHERE){
        return hittable_sphere_random((sphere*)o->obj, orig);
    }else if (o->t == HITTABLE_RECT){
        return hittable_rect_random((rect*)o->obj, orig);

    ///Wrappers
    }else if (o->t == HITTABLE_FLIP_FACE){
        return hittable_random(((flip_face*)o->obj)->obj, orig);

    }else{
        return (vec3){1,0,0};
    }
}

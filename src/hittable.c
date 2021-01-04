#include "hittable.h"



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
 * Hittable Objects allocators
 *
 * */


///Create a new sphere object
hittable* hittable_sphere_new(struct hittable_list* world, point3 center, double r, struct material* mat){
    //Create the sphere
    sphere* s = malloc(sizeof(sphere));
    s->center = center;
    s->r = r;
    s->mat = mat;

    //Create the hittable generic
    hittable* o = malloc(sizeof(hittable));
    o->obj = s;
    o->t = HITTABLE_SPHERE;

    //Add it to our world
    hittable_list_add(world, o);
    return o;
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

    //Create the hittable generic
    hittable* o = malloc(sizeof(hittable));
    o->obj = s;
    o->t = HITTABLE_MOVING_SPHERE;

    //Add it to our world
    hittable_list_add(world, o);
    return o;
}

point3 moving_sphere_center(moving_sphere* s, double t){
    double k = ((t - s->time0) / (s->time1 - s->time0));
    return vec3c_sum(vec3c_mul_k(vec3_sub(&s->center1, &s->center0), k), s->center0);
}




hittable* hittable_xy_rect_new(struct hittable_list* world, double x0, double x1, double y0, double y1, double k, struct material* mat){
    hittable* o = malloc(sizeof(hittable));
    o->t = HITTABLE_XYRECT;
    xy_rect* r = malloc(sizeof(xy_rect));
    r->x0 = x0; r->x1 = x1;
    r->y0 = y0; r->y1 = y1;
    r->k = k;
    r->mat = mat;
    o->obj = r;
    //Add it to our world
    hittable_list_add(world, o);
    return o;
}
hittable* hittable_xz_rect_new(struct hittable_list* world, double x0, double x1, double z0, double z1, double k, struct material* mat){
    hittable* o = malloc(sizeof(hittable));
    o->t = HITTABLE_XZRECT;
    xz_rect* r = malloc(sizeof(xz_rect));
    r->x0 = x0; r->x1 = x1;
    r->z0 = z0; r->z1 = z1;
    r->k = k;
    r->mat = mat;
    o->obj = r;
    //Add it to our world
    hittable_list_add(world, o);
    return o;
}
hittable* hittable_yz_rect_new(struct hittable_list* world, double y0, double y1, double z0, double z1, double k, struct material* mat){
    hittable* o = malloc(sizeof(hittable));
    o->t = HITTABLE_YZRECT;
    yz_rect* r = malloc(sizeof(yz_rect));
    r->y0 = y0; r->y1 = y1;
    r->z0 = z0; r->z1 = z1;
    r->k = k;
    r->mat = mat;
    o->obj = r;
    //Add it to our world
    hittable_list_add(world, o);
    return o;
}



hittable* hittable_box_new(struct hittable_list* world, point3 min, point3 max, struct material* mat){
    hittable* o = malloc(sizeof(hittable));
    o->t = HITTABLE_BOX;
    box* b = malloc(sizeof(box));
    b->box_max = max;
    b->box_min = min;

    //Make all 6 sides
    b->sides = (struct hittable_list*)hittable_list_new(6);
    hittable_xy_rect_new(b->sides, min.x, max.x, min.y, max.y, max.z, mat);
    hittable_xy_rect_new(b->sides, min.x, max.x, min.y, max.y, min.z, mat);

    hittable_xz_rect_new(b->sides, min.x, max.x, min.z, max.z, max.y, mat);
    hittable_xz_rect_new(b->sides, min.x, max.x, min.z, max.z, min.y, mat);

    hittable_yz_rect_new(b->sides, min.y, max.y, min.z, max.z, max.x, mat);
    hittable_yz_rect_new(b->sides, min.y, max.y, min.z, max.z, min.x, mat);

    o->obj = b;

    //Add it to our world
    hittable_list_add(world, o);
    return o;
}











int hittable_bounding_box(hittable* a, double t0, double t1, aabb* output_box);
inline static int box_compare(hittable* a, hittable* b, int axis){
    aabb box_a, box_b;
    if(!hittable_bounding_box(a, 0, 0, &box_a) || !hittable_bounding_box(b, 0, 0, &box_b)){
        printf("Failed to get bounding box in bvh_node comparator.\n"); fflush(stdout);
        return 0;
    }

    if (axis==0){       return box_a.minimum.x < box_b.minimum.x;
    }else if (axis==1){ return box_a.minimum.y < box_b.minimum.y;
    }else{              return box_a.minimum.z < box_b.minimum.z;}
}

int box_x_compare(const void* a, const void* b){return box_compare(*((hittable**)a), *((hittable**)b), 0);}
int box_y_compare(const void* a, const void* b){return box_compare(*((hittable**)a), *((hittable**)b), 1);}
int box_z_compare(const void* a, const void* b){return box_compare(*((hittable**)a), *((hittable**)b), 2);}




hittable* bvh_node_constructor(hittable** src_objects, int src_objects_count, int start, int end, double t0, double t1){
    hittable* h = malloc(sizeof(hittable));
    h->t = HITTABLE_BVH_NODE;

    bvh_node* node = malloc(sizeof(bvh_node));

    hittable** objects = malloc(sizeof(hittable*)*src_objects_count);
    memcpy(objects, src_objects, sizeof(hittable*)*src_objects_count);

    //printf("TEST: %p\n", objects[110]->obj);


    int raxis = random_double();

    int(*comparator)(const void*, const void*);
    if (raxis < 0.33){comparator = &box_x_compare;
    }else if (raxis < 0.66){comparator = &box_y_compare;
    }else{comparator = &box_z_compare;}

    int object_span = end-start;

    if(object_span == 1){
        node->left = objects[start];
        node->right = node->left;

    }else if (object_span == 2){
        if(comparator(&objects[start], &objects[start+1])){
            node->left = objects[start];
            node->right = objects[start+1];
        }else{
            node->left = objects[start+1];
            node->right = objects[start];
        }

    }else{
        qsort(objects, src_objects_count, sizeof(objects[0]), comparator);

        int mid = start + object_span/2;
        node->left  = bvh_node_constructor(src_objects, src_objects_count, start, mid, t0, t1);
        node->right = bvh_node_constructor(src_objects, src_objects_count, mid, end, t0, t1);
    }

    aabb box_left, box_right;
    if(!hittable_bounding_box(node->left, t0, t1, &box_left) || !hittable_bounding_box(node->right, t0, t1, &box_right)){
        printf("Failed to get bounding box in bvh_node constructor.\n"); fflush(stdout);
        return NULL;
    }
    node->box = surrounding_box(box_left, box_right);

    free(objects);

    h->obj = node;
    return h;
}


hittable* bvh_node_init(struct hittable_list* world, struct hittable_list* objects, double t0, double t1){
    void* objs = hittable_list_objects(objects);
    int index = hittable_list_index(objects);
    hittable* o = bvh_node_constructor(objs, index, 0, index, t0, t1);

    hittable_list_add(world, o);

    return o;
}







hittable* hittable_translate_init(struct hittable_list* world, hittable* obj, vec3 offset){
    hittable* o = malloc(sizeof(hittable));
    o->t = HITTABLE_TRANSLATE;
    translate* translated = malloc(sizeof(translate));
    translated->obj = obj;
    translated->offset = offset;
    o->obj = translated;
    hittable_list_add(world, o);
    return o;
}






hittable* hittable_rotate_y_init(struct hittable_list* world, hittable* obj, double angle){
    hittable* o = malloc(sizeof(hittable));
    o->t = HITTABLE_ROTATEY;
    rotate_y* rotated = malloc(sizeof(rotate_y));
    rotated->obj = obj;

    double radians = deg2rad(angle);
    rotated->sin_theta = sin(radians);
    rotated->cos_theta = cos(radians);
    rotated->hasbox = hittable_bounding_box(obj, 0, 1, &rotated->bbox);

    point3 min = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
    point3 max = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            for(int k=0; k<2; k++){
                double x = i * rotated->bbox.maximum.x + (1-i) * rotated->bbox.minimum.x;
                double y = j * rotated->bbox.maximum.y + (1-j) * rotated->bbox.minimum.y;
                double z = k * rotated->bbox.maximum.z + (1-k) * rotated->bbox.minimum.z;

                double newx = rotated->cos_theta * x + rotated->sin_theta * z;
                double newz = -rotated->sin_theta * x + rotated->cos_theta * z;

                vec3 tester = {newx, y, newz};
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

    o->obj = rotated;
    hittable_list_add(world, o);
    return o;
}









hittable* hittable_constant_medium_init(struct hittable_list* world, hittable* b, double d, texture* a){
    hittable* o = malloc(sizeof(hittable));
    o->t = HITTABLE_CONSTANT_MEDIUM;
    constant_medium* m = malloc(sizeof(constant_medium));
    m->boundary = b;
    m->neg_inv_density = -1/d;
    m->phase_function = (struct material*)material_isotropic_new(a);
    o->obj = m;
    hittable_list_add(world, o);
    return o;}
hittable* hittable_constant_medium_init_c(struct hittable_list* world, hittable* b, double d, color c){
    texture* t = texture_solid_color_init(c);
    return hittable_constant_medium_init(world, b, d, t);
}


















/**
 *
 * Hittable Objects deallocators
 *
 * */

///Free hittable object
void hittable_free(hittable* o){
    //printf("Freeing hittable %p (%d)\n", o, o->t);
    if(o->t == HITTABLE_SPHERE){
        sphere* s = o->obj;
        free(s);

    }else if(o->t == HITTABLE_MOVING_SPHERE){
        moving_sphere* s = o->obj;
        free(s);

    }else if(o->t == HITTABLE_BVH_NODE){

        bvh_node* n = o->obj;
        //printf("Freeing bvh leafs %p %p\n", n->left, n->right);

        if (n->left == n->right){
            hittable_free(n->left);
        }else{
            hittable_free(n->left);
            hittable_free(n->right);
        }
        free(n);

    }else if(o->t == HITTABLE_XYRECT){
        xy_rect* r = o->obj; free(r);
    }else if(o->t == HITTABLE_XZRECT){
        xz_rect* r = o->obj; free(r);
    }else if(o->t == HITTABLE_YZRECT){
        yz_rect* r = o->obj; free(r);

    }else if(o->t == HITTABLE_BOX){
        box* r = o->obj;
        free(r);

    }else if(o->t == HITTABLE_TRANSLATE){
        translate* r = o->obj;
        hittable_free(r->obj);
        free(r);

    }else if(o->t == HITTABLE_ROTATEY){
        rotate_y* r = o->obj;
        hittable_free(r->obj);
        free(r);

    }else if(o->t == HITTABLE_CONSTANT_MEDIUM){
        constant_medium* m = o->obj;
        hittable_free(m->boundary);
        free(m);

    }else{printf("!! Not implemented hittable_free for type %d\n", o->t);}
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
int hittable_hit(hittable* o, ray* r, double tmin, double tmax, hit_record* rec);

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



int bvh_node_hit(bvh_node* node, ray* r, double tmin, double tmax, hit_record* rec){
    if(aabb_hit(&node->box, r, tmin, tmax) == 0){return 0;}
    int hit_left = hittable_hit(node->left, r, tmin, tmax, rec);
    int hit_right = hittable_hit(node->right, r, tmin, hit_left ? rec->t : tmax, rec);
    return hit_left || hit_right;
}



int xy_rect_hit(xy_rect* rect, ray* r, double tmin, double tmax, hit_record* rec){
    double t = (rect->k - r->orig.z) / r->dir.z;
    if (t < tmin || t > tmax){return 0;}
    double x = r->orig.x + t * r->dir.x;
    double y = r->orig.y + t * r->dir.y;
    if (x < rect->x0 || x> rect->x1 || y < rect->y0 || y > rect->y1){return 0;}
    rec->u = (x - rect->x0)/(rect->x1 - rect->x0);
    rec->v = (y - rect->y0)/(rect->y1 - rect->y0);
    rec->t = t;
    vec3 outward_normal = {0,0,1};
    hit_record_set_facenormal(rec, r, &outward_normal);
    rec->mat = rect->mat;
    rec->p = ray_at(r, t);
    return 1;
}
int xz_rect_hit(xz_rect* rect, ray* r, double tmin, double tmax, hit_record* rec){
    double t = (rect->k - r->orig.y) / r->dir.y;
    if (t < tmin || t > tmax){return 0;}
    double x = r->orig.x + t * r->dir.x;
    double z = r->orig.z + t * r->dir.z;
    if (x < rect->x0 || x> rect->x1 || z < rect->z0 || z > rect->z1){return 0;}
    rec->u = (x - rect->x0)/(rect->x1 - rect->x0);
    rec->v = (z - rect->z0)/(rect->z1 - rect->z0);
    rec->t = t;
    vec3 outward_normal = {0,1,0};
    hit_record_set_facenormal(rec, r, &outward_normal);
    rec->mat = rect->mat;
    rec->p = ray_at(r, t);
    return 1;
}
int yz_rect_hit(yz_rect* rect, ray* r, double tmin, double tmax, hit_record* rec){
    double t = (rect->k - r->orig.x) / r->dir.x;
    if (t < tmin || t > tmax){return 0;}
    double y = r->orig.y + t * r->dir.y;
    double z = r->orig.z + t * r->dir.z;
    if (z < rect->z0 || z> rect->z1 || y < rect->y0 || y > rect->y1){return 0;}
    rec->u = (y - rect->y0)/(rect->y1 - rect->y0);
    rec->v = (z - rect->z0)/(rect->z1 - rect->z0);
    rec->t = t;
    vec3 outward_normal = {1,0,0};
    hit_record_set_facenormal(rec, r, &outward_normal);
    rec->mat = rect->mat;
    rec->p = ray_at(r, t);
    return 1;
}


int box_hit(box* b, ray* r, double tmin, double tmax, hit_record* rec){
    return hittable_list_hit(b->sides, r, tmin, tmax, rec);
}

int translate_hit(translate* b, ray* r, double tmin, double tmax, hit_record* rec){
    ray moved_ray = {vec3c_sub(r->orig, b->offset), r->dir, r->time};
    if (!hittable_hit(b->obj, &moved_ray, tmin, tmax, rec))
        return 0;

    rec->p = vec3c_sum(rec->p, b->offset);
    hit_record_set_facenormal(rec, &moved_ray, &rec->normal);
    return 1;
}


int rotate_y_hit(rotate_y* ob, ray* r, double tmin, double tmax, hit_record* rec){
    point3 origin = r->orig;
    vec3 direction = r->dir;

    origin.x = ob->cos_theta * r->orig.x  -  ob->sin_theta * r->orig.z;
    origin.z = ob->sin_theta * r->orig.x  +  ob->cos_theta * r->orig.z;

    direction.x = ob->cos_theta * r->dir.x  -  ob->sin_theta * r->dir.z;
    direction.z = ob->sin_theta * r->dir.x  +  ob->cos_theta * r->dir.z;

    ray rotated_ray = {origin, direction, r->time};
    if (!hittable_hit(ob->obj, &rotated_ray, tmin, tmax, rec))
        return 0;

    vec3 p = rec->p;
    vec3 normal = rec->normal;

    p.x = ob->cos_theta * rec->p.x  +  ob->sin_theta * rec->p.z;
    p.z = -ob->sin_theta * rec->p.x  +  ob->cos_theta * rec->p.z;

    normal.x = ob->cos_theta * rec->normal.x  +  ob->sin_theta * rec->normal.z;
    normal.z = -ob->sin_theta * rec->normal.x  +  ob->cos_theta * rec->normal.z;

    rec->p = p;
    hit_record_set_facenormal(rec, &rotated_ray, &normal);

    return 1;
}



int constant_medium_hit(constant_medium* m, ray* r, double tmin, double tmax, hit_record* rec){
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

    rec->t = rec1.t + hit_distance / ray_length;
    rec->p = ray_at(r, rec->t);

    rec->normal = (vec3){1,0,0};
    rec->front_face = 1;
    rec->mat = m->phase_function;
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
    }else if (o->t == HITTABLE_XYRECT){
        return xy_rect_hit((xy_rect*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_XZRECT){
        return xz_rect_hit((xz_rect*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_YZRECT){
        return yz_rect_hit((yz_rect*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_BOX){
        return box_hit((box*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_TRANSLATE){
        return translate_hit((translate*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_ROTATEY){
        return rotate_y_hit((rotate_y*)o->obj, r, tmin, tmax, rec);
    }else if (o->t == HITTABLE_CONSTANT_MEDIUM){
        return constant_medium_hit((constant_medium*)o->obj, r, tmin, tmax, rec);
    }else{printf("!! Not implemented hittable_hit for type %d\n", o->t);return 0;}
}





/**
 *
 * Hittable Objects Bounding box
 *
 * */


int sphere_bounding_box(sphere* s, aabb* output_box){
    output_box->minimum = vec3c_sub(s->center, (vec3){s->r,s->r,s->r});
    output_box->maximum = vec3c_sum(s->center, (vec3){s->r,s->r,s->r});
    return 1;
}


int moving_sphere_bounding_box(moving_sphere* s, double t0, double t1, aabb* output_box){
    vec3 area = {s->r, s->r, s->r};
    aabb box0 = {vec3c_sub(moving_sphere_center(s, t0), area), vec3c_sum(moving_sphere_center(s, t0), area)};
    aabb box1 = {vec3c_sub(moving_sphere_center(s, t1), area), vec3c_sum(moving_sphere_center(s, t1), area)};
    aabb res = surrounding_box(box0, box1);
    output_box->minimum = res.minimum;
    output_box->maximum = res.maximum;
    return 1;
}



int bvh_node_bounding_box(bvh_node* node, aabb* output_box){
    *output_box = node->box;
    return 1;
}


int xy_rect_bounding_box(xy_rect* rect, aabb* output_box){
    point3 a = {rect->x0, rect->y0, rect->k-0.0001};
    point3 b = {rect->x1, rect->y1, rect->k+0.0001};
    *output_box = (aabb){a, b};
    return 1;}
int xz_rect_bounding_box(xz_rect* rect, aabb* output_box){
    point3 a = {rect->x0, rect->k-0.0001, rect->z0};
    point3 b = {rect->x1, rect->k+0.0001, rect->z1};
    *output_box = (aabb){a, b};
    return 1;}
int yz_rect_bounding_box(yz_rect* rect, aabb* output_box){
    point3 a = {rect->k-0.0001, rect->y0, rect->z0};
    point3 b = {rect->k+0.0001, rect->y1, rect->z1};
    *output_box = (aabb){a, b};
    return 1;}


int box_bounding_box(box* b, aabb* output_box){
    *output_box = (aabb){b->box_min, b->box_max};
    return 1;
}

int translate_bounding_box(translate* t, double t0, double t1, aabb* output_box){
    if (!hittable_bounding_box(t->obj, t0, t1, output_box)) return 0;
    *output_box = (aabb){vec3c_sum(output_box->minimum, t->offset), vec3c_sum(output_box->maximum, t->offset)};
    return 1;
}

int rotate_y_bounding_box(rotate_y* ob, aabb* output_box){
    *output_box = ob->bbox;
    return ob->hasbox;
}

int constant_medium_bounding_box(constant_medium* m, double t0, double t1, aabb* output_box){
    return hittable_bounding_box(m->boundary, t0, t1, output_box);
}



///Generic function to test for ray hit
int hittable_bounding_box(hittable* o, double t0, double t1, aabb* output_box){
    if (o->t == HITTABLE_SPHERE){
        return sphere_bounding_box((sphere*)o->obj, output_box);
    }else if (o->t == HITTABLE_MOVING_SPHERE){
        return moving_sphere_bounding_box((moving_sphere*)o->obj, t0, t1, output_box);
    }else if (o->t == HITTABLE_BVH_NODE){
        return bvh_node_bounding_box((bvh_node*)o->obj, output_box);
    }else if (o->t == HITTABLE_XYRECT){
        return xy_rect_bounding_box((xy_rect*)o->obj, output_box);
    }else if (o->t == HITTABLE_XZRECT){
        return xz_rect_bounding_box((xz_rect*)o->obj, output_box);
    }else if (o->t == HITTABLE_YZRECT){
        return yz_rect_bounding_box((yz_rect*)o->obj, output_box);
    }else if (o->t == HITTABLE_BOX){
        return box_bounding_box((box*)o->obj, output_box);
    }else if (o->t == HITTABLE_TRANSLATE){
        return translate_bounding_box((translate*)o->obj, t0, t1, output_box);
    }else if (o->t == HITTABLE_ROTATEY){
        return rotate_y_bounding_box((rotate_y*)o->obj, output_box);
    }else if (o->t == HITTABLE_CONSTANT_MEDIUM){
        return constant_medium_bounding_box((constant_medium*)o->obj, t0, t1, output_box);
    }else{printf("!! Not implemented hittable_hit for type %d\n", o->t);return 0;}
}

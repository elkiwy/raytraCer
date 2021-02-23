/**
 * Utils
 **/
#define PI 3.14159265359

typedef struct Inputs {
   __global float16* objects;
   __global float16* wrapped_objects;
   __global float16* lights;
   __global float16* materials;
   __global float16* textures;
   uint objs_count;
   uint wrapped_objs_count;
   uint lights_count;
   uint mats_count;
   uint texs_count;

    uint allocated_objects;
    uint allocated_lights;

   __global float3* rus;
   int rus_count;
   int* seed;

   __global float3* perlin_ranvec;
   __global int3* perlin_perms;
   int perlin_point_count;

   int GLOBAL_ID;
   int ITERATION;
   int SAMPLE;
} Inputs;

float3 rus_to_rhemi(float3 rus, float3 normal);
int rand(int* seed);
float rand_float(int* seed);
float rand_float_scaled(int* seed, float min, float max);
float3 random_in_unit_sphere(const Inputs* i);
float3 random_in_unit_vector(const Inputs* i);
float3 random_in_hemisphere(const Inputs* i, float3 normal);




//Random in unit sphere to random in hemisphere
float3 rus_to_rhemi(float3 rus, float3 normal){
   if(dot(rus, normal)>0.0){
      return rus;
   }else{
      return -rus;
   }
}


///Generate random number between 0 and 2147483647
int rand(int* seed) {// 1 <= *seed < m
    int const a = 16807; //ie 7**5
    int const m = 2147483647; //ie 2**31-1
    int val = (long)(*(seed) * a) % m;
    //*seed = abs(val);

    if (val<0){*seed = -1*val; //TODO: optimize this, need to force positive on AMD GPU
    }else{*seed = val;}
    return(*seed);
}

///Generate random float number between 0-1
float rand_float(int* seed){
   float v = rand(seed);
   return (float)(v / 2147483647.0f);
}

float rand_float_scaled(int* seed, float min, float max){
   float v = rand_float(seed);
   return (v*(max-min))+min;
}

float3 random_in_unit_sphere(const Inputs* i){
   return i->rus[rand(i->seed) % i->rus_count];
}

float3 random_in_unit_vector(const Inputs* i){
   return normalize(random_in_unit_sphere(i));
}

float3 random_in_hemisphere(const Inputs* i, float3 normal){
    float3 rus = i->rus[rand(i->seed) % i->rus_count];
    if(dot(rus, normal)>0.0){
        return rus;
    }else{
        return -rus;
    }
}




typedef struct onb{
    float3 u;
    float3 v;
    float3 w;
}onb;

void onb_init(onb* o, const float3 n);
void onb_init(onb* o, const float3 n){
    o->w = normalize(n);
    float3 a = (fabs(o->w[0]) > 0.9) ? (float3){0,1,0} : (float3){1,0,0};
    o->v = normalize(cross(o->w, a));
    o->u = cross(o->w, o->v);
}

float3 onb_local(const onb* o, const float3 a);
float3 onb_local(const onb* o, const float3 a){
    return a[0]*o->u + a[1]*o->v + a[2]*o->w;
}




void print16(__constant char* s, float16 v);
void print16(__constant char* s, float16 v){
   printf("%s:\n", s);
   for (int i=0;i<16;i++){
      printf("  v[%d]: %f\n", i, v[i]);
   }
}


void print3(__constant char* s, float3 v){
   printf("%s:\n", s);
   for (int i=0;i<3;i++){
      printf("  v[%d]: %f\n", i, v[i]);
   }
}










/**
 * Rays
 **/

typedef struct ray {
   float3 orig;
   float3 dir;
} ray;

float3 ray_at(const ray* r, float t);


float3 ray_at(const ray* r, float t) {
   return r->orig + t * r->dir;
}


void printray(__constant char* s, ray r){
    printf("%s: \n %f %f %f  ->  %f %f %f   control:%f\n",
           s,
           r.orig[0], r.orig[1], r.orig[2],
           r.dir[0], r.dir[1], r.dir[2],
           123.0f);
}

















/**
 * Hittable
 **/

//struct material;

typedef struct hit_record {
   float3 p;
   float3 normal;
   float t;
   int front_face;
   int mat_ptr;
   float u, v;
} hit_record;

void set_face_normal(hit_record* rec, const ray* r, float3 outward_normal);

void set_face_normal(hit_record* rec, const ray* r, float3 outward_normal) {
   rec->front_face = dot(r->dir, outward_normal) < 0;
   rec->normal = rec->front_face ? outward_normal :-outward_normal;
}

//typedef struct hittable {
//   uint type;
//   float3 p1;
//   float3 p2;
//   float radius;
//   struct material* mat;
//} hittable;

#define OBJK_TYPE 0
#define OBJK_X0 1
#define OBJK_Y0 2
#define OBJK_Z0 3
#define OBJK_X1 4
#define OBJK_Y1 5
#define OBJK_Z1 6
#define OBJK_RAD  7
#define OBJK_AXIS 7
#define OBJK_MAT  8
#define OBJK_RECT_K  9
#define OBJK_DENSITY 9
#define OBJK_ROTATION_ANGLE  9

#define OBJK_OBJPTR1  10
#define OBJK_OBJPTR2  11
#define OBJK_OBJPTR3  12
#define OBJK_OBJPTR4  13
#define OBJK_OBJPTR5  14
#define OBJK_OBJPTR6  15

float length_squared(float3 v);
bool hittable_hit(const Inputs* inputs, const float16 obj, const ray* r, double t_min, double t_max, hit_record* rec);
bool worldHit(const Inputs* inputs, const __global float16* objs, const int objs_count, const ray* r, hit_record *rec);

float length_squared(float3 v){return v.x*v.x + v.y*v.y + v.z*v.z;}




void hittable_get_uv(const float16 obj, float3 p, float* u, float* v);
void hittable_get_uv(const float16 obj, float3 p, float* u, float* v){
   if(obj[OBJK_TYPE]==0){

      float theta = acos(-p[1]);
      float phi = atan2(-p[2], p[0]) + PI;
      *u = phi / (2*PI);
      *v = theta / PI;

   }else{
      //TODO
   }
}


#define XY 0
#define YZ 1
#define XZ 2


///Rotate a point in space, useful on hittable_rotate
void rotatePoint(float3 original, const float sin_theta, const float cos_theta, const int axis, float3* output);
void rotatePoint(float3 original, const float sin_theta, const float cos_theta, const int axis, float3* output){
    if (axis == 0){
        output->x =  cos_theta * original.x  -  sin_theta * original.y;
        output->y =  sin_theta * original.x  +  cos_theta * original.y;
    }else if (axis == 1){
        output->x =  cos_theta * original.x  +  sin_theta * original.z;
        output->z = -sin_theta * original.x  +  cos_theta * original.z;
    }else{
        output->y =  cos_theta * original.y  -  sin_theta * original.z;
        output->z =  sin_theta * original.y  +  cos_theta * original.z;
    }
}


bool hittable_hit(const Inputs* inputs, const float16 obj_to_test, const ray* r, double t_min, double t_max, hit_record* rec) {

   float16 objects[16];
   int object_to_test = 1;
   objects[0] = obj_to_test;
   bool hit_anything = false;
   float closestSoFar = t_max;

   ray current_ray = *r;
   float3 ray_offset = {0,0,0};

   float rotation_sin_theta = 0;
   float rotation_cos_theta = 1;
   float rotation_axis = 0;
   int flipping_face = 0;

   for(int i=0;i<object_to_test;i++){

      float16 obj = objects[i];

      //if (i>0){printf("Testing %f %d\n", obj[OBJK_TYPE], i);}

      //Sphere
      if (obj[OBJK_TYPE]==0){
         float3 center = {obj[OBJK_X0],obj[OBJK_Y0],obj[OBJK_Z0]};
         float3 oc = current_ray.orig - center;
         float a = length_squared(current_ray.dir);
         float half_b = dot(oc, current_ray.dir);
         float c = length_squared(oc) - obj[OBJK_RAD]*obj[OBJK_RAD];
         float discriminant = half_b*half_b - a*c;
         if (discriminant < 0) {continue;}

         float sqrtd = sqrt(discriminant);
         float root  = (-half_b - sqrtd) / a;
         if (root < t_min || closestSoFar < root) {
            root = (-half_b + sqrtd) / a;
            if (root < t_min || closestSoFar < root){continue;}
         }

         rec->t = root;

         rec->p = ray_at(&current_ray, rec->t);
         rec->p += ray_offset;

         float3 outward_normal = (rec->p - center) / obj[OBJK_RAD];
         ray offsetted_ray = {current_ray.orig - ray_offset, current_ray.dir};
         set_face_normal(rec, &offsetted_ray, outward_normal);

         hittable_get_uv(obj, outward_normal, &rec->u, &rec->v);

         rec->mat_ptr = (int)round(obj[OBJK_MAT]);

         hit_anything = true;
         closestSoFar = rec->t;

      //Rect
      }else if (obj[OBJK_TYPE] == 1){

          //Check for plane ray intersection
          float t;
          if(obj[OBJK_AXIS] == XY){       t = (obj[OBJK_RECT_K] - current_ray.orig[2]) / current_ray.dir[2];
          }else if(obj[OBJK_AXIS] == XZ){ t = (obj[OBJK_RECT_K] - current_ray.orig[1]) / current_ray.dir[1];
          }else{                          t = (obj[OBJK_RECT_K] - current_ray.orig[0]) / current_ray.dir[0];}
          if (t < t_min || t > closestSoFar){
              continue;
          }

          //Find if the point of intersection is in the boundary
          float x = current_ray.orig[0] + t * current_ray.dir[0]; float y = current_ray.orig[1] + t * current_ray.dir[1]; float z = current_ray.orig[2] + t * current_ray.dir[2];
          if ((obj[OBJK_AXIS] == XZ && (x < obj[OBJK_X0] || x> obj[OBJK_X1] || z < obj[OBJK_Z0] || z > obj[OBJK_Z1]))
              || (obj[OBJK_AXIS] == XY && (x < obj[OBJK_X0] || x> obj[OBJK_X1] || y < obj[OBJK_Y0] || y > obj[OBJK_Y1]))
              || (obj[OBJK_AXIS] == YZ && (z < obj[OBJK_Z0] || z> obj[OBJK_Z1] || y < obj[OBJK_Y0] || y > obj[OBJK_Y1]))) {
              continue;
          }

          //Register the hit
          rec->t = t;
          rec->p = ray_at(&current_ray, rec->t);
          //Revert rotation on the normal
          float3 p = rec->p;
          float3 normal = rec->normal;
          rotatePoint(rec->p, rotation_sin_theta, rotation_cos_theta, rotation_axis, &p);
          rotatePoint(rec->normal, rotation_sin_theta, rotation_cos_theta, rotation_axis, &normal);
          rec->p = p;
          //Revert translation
          rec->p += ray_offset;


          float3 outward_normal = {obj[OBJK_AXIS]==YZ?1:0,obj[OBJK_AXIS]==XZ?1:0,obj[OBJK_AXIS]==XY?1:0};
          ray offsetted_ray = {current_ray.orig - ray_offset, current_ray.dir};
          set_face_normal(rec, &offsetted_ray, outward_normal);

          //Apply eventual flipping face
          rec->front_face = flipping_face ^ rec->front_face;

          rec->u = obj[OBJK_AXIS]==YZ ? (y - obj[OBJK_Y0])/(obj[OBJK_Y1] - obj[OBJK_Y0]) : (x - obj[OBJK_X0])/(obj[OBJK_X1] - obj[OBJK_X0]);
          rec->v = obj[OBJK_AXIS]==XY ? (y - obj[OBJK_Y0])/(obj[OBJK_Y1] - obj[OBJK_Y0]) : (z - obj[OBJK_Z0])/(obj[OBJK_Z1] - obj[OBJK_Z0]);

          rec->mat_ptr = (int)round(obj[OBJK_MAT]);

          hit_anything = true;
          closestSoFar = rec->t;


      //Box
      }else if (obj[OBJK_TYPE] == 2){


         object_to_test += 6;
         objects[i+1] = inputs->wrapped_objects[(int)round(obj[OBJK_OBJPTR1])];
         objects[i+2] = inputs->wrapped_objects[(int)round(obj[OBJK_OBJPTR2])];
         objects[i+3] = inputs->wrapped_objects[(int)round(obj[OBJK_OBJPTR3])];
         objects[i+4] = inputs->wrapped_objects[(int)round(obj[OBJK_OBJPTR4])];
         objects[i+5] = inputs->wrapped_objects[(int)round(obj[OBJK_OBJPTR5])];
         objects[i+6] = inputs->wrapped_objects[(int)round(obj[OBJK_OBJPTR6])];

      //Rotated
      }else if (obj[OBJK_TYPE] == 3){

          float angle = obj[OBJK_ROTATION_ANGLE] * PI/180.0f;
          rotation_sin_theta = sin(angle);
          rotation_cos_theta = cos(angle);
          rotation_axis = (int)round(obj[OBJK_AXIS]);

          //Rotate the incoming ray
          float3 rot_origin = current_ray.orig;
          float3 rot_direction = current_ray.dir;
          rotatePoint(current_ray.orig, -rotation_sin_theta, rotation_cos_theta, rotation_axis, &rot_origin);
          rotatePoint(current_ray.dir, -rotation_sin_theta, rotation_cos_theta, rotation_axis, &rot_direction);

          //Test on the wrapped object
          current_ray = (ray){rot_origin, rot_direction};
          object_to_test += 1;
          objects[i+1] = inputs->wrapped_objects[(int)round(obj[OBJK_OBJPTR1])];



      //Translated
      }else if (obj[OBJK_TYPE] == 4){

          ray_offset = (float3){obj[OBJK_X0], obj[OBJK_Y0], obj[OBJK_Z0]};
          current_ray = (ray){current_ray.orig - ray_offset, current_ray.dir};

          object_to_test += 1;
          objects[i+1] = inputs->wrapped_objects[(int)round(obj[OBJK_OBJPTR1])];

      //Constant Medium
      }else if (obj[OBJK_TYPE] == 5){
          //SUPPORT ONLY SPHERES
          float neg_inv_density = obj[OBJK_DENSITY];

          //Get the two intersection from the wrapped sphere
          float3 center = {obj[OBJK_X0],obj[OBJK_Y0],obj[OBJK_Z0]};
          float3 oc = current_ray.orig - center;
          float a = length_squared(current_ray.dir);
          float half_b = dot(oc, current_ray.dir);
          float c = length_squared(oc) - obj[OBJK_RAD]*obj[OBJK_RAD];
          float discriminant = half_b*half_b - a*c;
          if (discriminant < 0) {continue;}
          float sqrtd = sqrt(discriminant);
          float root1 = (-half_b - sqrtd) / a;
          float root2 = (-half_b + sqrtd) / a;
          if (root1 < t_min || closestSoFar < root1 || root2 < t_min || closestSoFar < root2){continue;}
          if (root1 < t_min) root1 = t_min;
          if (root2 < t_max) root2 = t_max;
          if (root1 >= root2){continue;}
          root1 = max(0.0f, root1);

          const float ray_length = sqrt(length_squared(current_ray.dir));
          const float distance_inside_sphere = (root2 - root1) * ray_length;
          const float hit_distance = neg_inv_density * log(rand_float(inputs->seed));
          if (hit_distance > distance_inside_sphere){continue;}
          rec->t = root1 + hit_distance / ray_length;
          rec->p = ray_at(&current_ray, rec->t);
          rec->normal = (float3){1,0,0};
          rec->front_face = 1;
          rec->mat_ptr = (int)round(obj[OBJK_MAT]);

          hit_anything = true;
          closestSoFar = rec->t;

      //Flip face
      }else if (obj[OBJK_TYPE] == 6){
          flipping_face = 1;

          object_to_test += 1;
          objects[i+1] = inputs->wrapped_objects[(int)round(obj[OBJK_OBJPTR1])];
      }
   }
   return hit_anything;
}





bool worldHit(const Inputs* inputs, const __global float16* objs, const int objs_count, const ray* r, hit_record *rec) {
   hit_record temp_rec;
   bool hitAnything = false;
   float closestSoFar = FLT_MAX;

   for (int i = 0; i < objs_count; i++) {
      if (hittable_hit(inputs, objs[i], r, 0.001, closestSoFar, &temp_rec)) {
         hitAnything = true;
         closestSoFar = temp_rec.t;
         *rec = temp_rec;
      }
   }

   return hitAnything;
}






double hittable_pdf_value(const Inputs* inputs, const int obj_ptr, const float3 o, const float3 v);
double hittable_pdf_value(const Inputs* inputs, const int obj_ptr, const float3 o, const float3 v){
    float16 obj = inputs->lights[obj_ptr];
    //Sphere
    if (obj[OBJK_TYPE] == 0){
        hit_record rec;
        ray r = {o, v};
        if(!hittable_hit(inputs, obj, &r, 0.0001, FLT_MAX, &rec)){return 0;}

        float3 center = {obj[OBJK_X0],obj[OBJK_Y0],obj[OBJK_Z0]};
        float cos_theta_max = sqrt(1.0f-obj[OBJK_RAD]*obj[OBJK_RAD]/length_squared(center-o));
        float solid_angle = 2.0f*PI*(1.0f-cos_theta_max);
        return 1 / solid_angle;

    //Rect
    }else if (obj[OBJK_TYPE] == 1){
        hit_record rec;
        ray r = {o, v};
        if(!hittable_hit(inputs, obj, &r, 0.0001, FLT_MAX, &rec)){return 0;}
        float distance_squared = rec.t * rec.t * length_squared(v);

        float area = 0;
        if(      obj[OBJK_AXIS] == XZ){area = (obj[OBJK_X1]-obj[OBJK_X0])*(obj[OBJK_Z1]-obj[OBJK_Z0]);
        }else if(obj[OBJK_AXIS] == XY){area = (obj[OBJK_X1]-obj[OBJK_X0])*(obj[OBJK_Y1]-obj[OBJK_Y0]);
        }else if(obj[OBJK_AXIS] == YZ){area = (obj[OBJK_Y1]-obj[OBJK_Y0])*(obj[OBJK_Z1]-obj[OBJK_Z0]);}
        float cosine = fabs(dot(v, rec.normal) / length(v));
        return distance_squared / (cosine * area);
    }

    return 0.0;
}

float3 random_to_sphere(const Inputs* inputs, float radius, float distance_squared);
float3 random_to_sphere(const Inputs* inputs, float radius, float distance_squared) {
    float r1 = rand_float(inputs->seed);
    float r2 = rand_float(inputs->seed);
    float z = 1.0f + r2*(sqrt(1.0f-radius*radius/distance_squared) - 1.0f);
    float phi = 2.0f*PI*r1;
    float x = cos(phi)*sqrt(1.0f-z*z);
    float y = sin(phi)*sqrt(1.0f-z*z);
    return (float3){x, y, z};
}

float3 hittable_random(const Inputs* inputs, const int obj_ptr, const float3 o);
float3 hittable_random(const Inputs* inputs, const int obj_ptr, const float3 o){
    float16 obj = inputs->lights[obj_ptr];

    //Sphere
    if (obj[OBJK_TYPE] == 0){
        float3 center = {obj[OBJK_X0],obj[OBJK_Y0],obj[OBJK_Z0]};
        float3 direction = center - o;
        float distance_squared = length_squared(direction);
        onb uvw;
        onb_init(&uvw, direction);
        return onb_local(&uvw, random_to_sphere(inputs, obj[OBJK_RAD], distance_squared));

    //Rect
    }if (obj[OBJK_TYPE] == 1){

        float3 random_point;
        if(obj[OBJK_AXIS] == XZ){
            random_point = (float3){rand_float_scaled(inputs->seed, obj[OBJK_X0], obj[OBJK_X1]), obj[OBJK_RECT_K], rand_float_scaled(inputs->seed, obj[OBJK_Z0], obj[OBJK_Z1])};
        }else if(obj[OBJK_AXIS] == XY){
            random_point = (float3){rand_float_scaled(inputs->seed, obj[OBJK_X0], obj[OBJK_X1]), rand_float_scaled(inputs->seed, obj[OBJK_Y0], obj[OBJK_Y1]), obj[OBJK_RECT_K]};
        }else{
            random_point = (float3){obj[OBJK_RECT_K], rand_float_scaled(inputs->seed, obj[OBJK_Y0], obj[OBJK_Y1]), rand_float_scaled(inputs->seed, obj[OBJK_Z0], obj[OBJK_Z1])};
        }
        return random_point - o;
    }
    return (float3){1,0,0};
}







/**
 * PDF
 **/

#define PDFT_COSINE 1
#define PDFT_OBJECT 2
#define PDFT_MIXTURE 3

typedef struct pdf{
    //Generic
    int type; const Inputs* inputs;
    //Cosine
    onb uvw;
    //Object
    int obj_ptr; float3 o;
    //Mixture
}pdf;




float pdf_value(const pdf* p, const float3 direction);
float pdf_value(const pdf* p, const float3 direction){
    //Cosine
    if (p->type == PDFT_COSINE){
        float cosine = dot(normalize(direction), p->uvw.w);
        return (cosine <= 0) ? 0 : cosine/PI;

    //Object
    }else if (p->type == PDFT_OBJECT){
        return hittable_pdf_value(p->inputs, p->obj_ptr, p->o, direction);

    //Mixture
    }else if (p->type == PDFT_MIXTURE){
        float cosine = dot(normalize(direction), p->uvw.w);
        float cos_val = (cosine <= 0) ? 0 : cosine/PI;
        float obj_val = hittable_pdf_value(p->inputs, p->obj_ptr, p->o, direction);
        return 0.5f * cos_val + 0.5f * obj_val;
    }
    return 0.0;
}



float3 random_cosine_direction(const Inputs* inputs);
float3 random_cosine_direction(const Inputs* inputs) {
    float r1 = rand_float(inputs->seed);
    float r2 = rand_float(inputs->seed);
    float z = sqrt(1-r2);
    float phi = 2*PI*r1;
    float x = cos(phi)*sqrt(r2);
    float y = sin(phi)*sqrt(r2);
    return (float3){x, y, z};
}
float3 pdf_generate(const pdf* p);
float3 pdf_generate(const pdf* p){
    //Cosine
    if (p->type == PDFT_COSINE){
        return onb_local(&p->uvw, random_cosine_direction(p->inputs));

    //Object
    }else if (p->type == PDFT_OBJECT){
        return hittable_random(p->inputs, p->obj_ptr, p->o);

    //Mixture
    }else if (p->type == PDFT_MIXTURE){
        if(rand_float(p->inputs->seed) < 0.5){
            //Object
            return hittable_random(p->inputs, p->obj_ptr, p->o);
        }else{
            //Cosine
            return onb_local(&p->uvw, random_cosine_direction(p->inputs));
        }
    }
    return (float3){1,0,0};
}
















/**
 * Textures
 **/


float noise(const Inputs* inp, float3 p);
float turb(const Inputs* inp, float3 p, int depth);

float noise(const Inputs* inp, float3 p) {
   float u = p[0] - floor(p[0]);
   float v = p[1] - floor(p[1]);
   float w = p[2] - floor(p[2]);
   const int i = (int)floor(p[0]);
   const int j = (int)floor(p[1]);
   const int k = (int)floor(p[2]);

   float3 c[2][2][2];
   for(int di=0;di<2;di++){
      for(int dj=0;dj<2;dj++){
         for(int dk=0;dk<2;dk++){
            c[di][dj][dk] = inp->perlin_ranvec[
               inp->perlin_perms[(i+di) & 255][0] ^
               inp->perlin_perms[(j+dj) & 255][1] ^
               inp->perlin_perms[(k+dk) & 255][2]
            ];
         }
      }
   }

   float uu = u*u*(3-2*u);
   float vv = v*v*(3-2*v);
   float ww = w*w*(3-2*w);
   float accum = 0.0;

   for (int i=0; i < 2; i++)
      for (int j=0; j < 2; j++)
         for (int k=0; k < 2; k++) {
            float3 weight_v = {u-i, v-j, w-k};
            accum += (i*uu + (1-i)*(1-uu))
               * (j*vv + (1-j)*(1-vv))
               * (k*ww + (1-k)*(1-ww))
               * dot(c[i][j][k], weight_v);
         }

   return accum;
}

float turb(const Inputs* inp, float3 p, int depth) {
   float accum = 0.0;
   float3 temp_p = p;
   float weight = 1.0;

   for (int i = 0; i < depth; i++) {
      accum += weight * noise(inp, temp_p);
      weight *= 0.5;
      temp_p *= 2;
   }

   return fabs(accum);
}

#define TEXK_TYPE 0
#define TEXK_COL0 1
#define TEXK_COL1 2
#define TEXK_COL2 3
#define TEXK_SCALE 4

float3 texture_value(const Inputs* inputs, const __global float16* tex_ptr, float u, float v, float3 p);
float3 texture_value(const Inputs* inputs, const __global float16* tex_ptr, float u, float v, float3 p){
   float16 tex = *tex_ptr;
   if (tex[TEXK_TYPE] == 0){
      return (float3){tex[TEXK_COL0], tex[TEXK_COL1], tex[TEXK_COL2]};
   }else if(tex[TEXK_TYPE] == 1){
      float v = 0.5 * (1.0f + sin(turb(inputs, p, 7)*10 + tex[TEXK_SCALE]*p[2]));
      return (float3){1,1,1} * v;
   }
   return (float3){0,0,0};
}











/**
 * Materials
 **/

typedef struct scatter_record {
   ray specular_ray;
   int is_specular;
   float3 attenuation;
   pdf scatter_pdf;
} scatter_record;


#define MATK_TYPE 0
#define MATK_ALB0 1
#define MATK_ALB1 2
#define MATK_ALB2 3
#define MATK_FUZZ 4
#define MATK_IR 5
#define MATK_TEX_PTR 6

bool near_zero(float3 v);
bool near_zero(float3 v){
   const float s = 1e-8;
   return (fabs(v[0])<s) && (fabs(v[1])<s) && (fabs(v[2])<s);
}

float3 reflect(float3 v, float3 n);
float3 reflect(float3 v, float3 n){
   return v - 2*dot(v, n)*n;
}

float3 refract(float3 uv, float3 n, float etai_over_etat);
float3 refract(float3 uv, float3 n, float etai_over_etat){
   float cos_theta = fmin((float)dot(-uv, n), 1.0f);
   float3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
   float3 r_out_parallel = -sqrt(fabs(1.0f - length_squared(r_out_perp))) * n;
   return r_out_perp + r_out_parallel;
}


float reflectance(float cosine, float ref_idx) ;
float reflectance(float cosine, float ref_idx) {
   // Use Schlick's approximation for reflectance.
   float r0 = (1.0f-ref_idx) / (1.0f+ref_idx);
   r0 = r0*r0;
   return r0 + (1.0f-r0)*pow((1.0f - cosine),5.0f);
}



float3 material_emitted(const Inputs* inputs, const int mat_ptr, const ray* r_in, const hit_record* rec,  double u, double v, float3 p);
float3 material_emitted(const Inputs* inputs, const int mat_ptr, const ray* r_in, const hit_record* rec,  double u, double v, float3 p){
    float16 mat = inputs->materials[mat_ptr];
    //Light
    if (mat[MATK_TYPE] == 3){
        if (rec->front_face == 1){
            return texture_value(inputs, &(inputs->textures[(int)round(mat[MATK_TEX_PTR])]), u, v, p);
        }
    }
    return (float3){0,0,0};
}



void pdf_init_cosine(pdf* p, const float3 normal, const Inputs* inputs);
void pdf_init_cosine(pdf* p, const float3 normal, const Inputs* inputs){
    p->type = PDFT_COSINE;
    p->inputs = inputs;
    onb_init(&p->uvw, normal);
}

void pdf_init_object(pdf* p, const int obj_ptr, const float3 orig, const Inputs* inputs);
void pdf_init_object(pdf* p, const int obj_ptr, const float3 orig, const Inputs* inputs){
    p->type = PDFT_OBJECT;
    p->o = orig;
    p->obj_ptr = obj_ptr;
    p->inputs = inputs;
}

void pdf_init_mixture(pdf* p, const float3 normal, const int obj_ptr, const float3 orig, const Inputs* inputs);
void pdf_init_mixture(pdf* p, const float3 normal, const int obj_ptr, const float3 orig, const Inputs* inputs){
    p->type = PDFT_MIXTURE;
    p->o = orig;
    p->obj_ptr = obj_ptr;
    p->inputs = inputs;
    onb_init(&p->uvw, normal);
}



bool material_scatter(const Inputs* inputs, const int mat_ptr, const ray* r_in, const hit_record* rec, scatter_record* srec) ;
bool material_scatter(const Inputs* inputs, const int mat_ptr, const ray* r_in, const hit_record* rec, scatter_record* srec) {
   float16 mat = inputs->materials[mat_ptr];

   //Lambertian
   if (mat[MATK_TYPE]==0){
       srec->is_specular = 0;
       srec->attenuation = texture_value(inputs, &(inputs->textures[(int)round(mat[MATK_TEX_PTR])]), rec->u, rec->v, rec->p);
       pdf_init_cosine(&srec->scatter_pdf, rec->normal, inputs);
       return true;

   //Metal
   }else if (mat[MATK_TYPE] == 1){
      float3 reflected = reflect(normalize(r_in->dir), rec->normal);
      srec->specular_ray = (ray){rec->p, reflected + (mat[MATK_FUZZ] * random_in_unit_sphere(inputs))};
      srec->attenuation = texture_value(inputs, &(inputs->textures[(int)round(mat[MATK_TEX_PTR])]), rec->u, rec->v, rec->p);
      srec->is_specular = 1;
      return true;

   //Dielectric
   }else if (mat[MATK_TYPE] == 2){

       srec->is_specular = 1;
       srec->attenuation = (float3){1.0,1.0,1.0};

       float refraction_ratio = rec->front_face ? (1.0f/mat[MATK_IR]) : mat[MATK_IR];
       float3 unit_direction = normalize(r_in->dir);
       float cos_theta = fmin(dot(-unit_direction, rec->normal), 1.0f);
       float sin_theta = sqrt(1.0f - cos_theta*cos_theta);
       bool cannot_refract = refraction_ratio * sin_theta > 1.0f;
       float3 direction;
       if (cannot_refract || reflectance(cos_theta, refraction_ratio) > rand_float(inputs->seed)){
           direction = reflect(unit_direction, rec->normal);
       }else{
           direction = refract(unit_direction, rec->normal, refraction_ratio);
       }
       srec->specular_ray = (ray){rec->p, direction};
       return true;

    //Light
   }else if (mat[MATK_TYPE] == 3){
      return false;

   //Isotropic
   }else if (mat[MATK_TYPE] == 4){
       srec->is_specular = 1;
       srec->specular_ray = (ray){rec->p, random_in_unit_sphere(inputs)};
       srec->attenuation = texture_value(inputs, &(inputs->textures[(int)round(mat[MATK_TEX_PTR])]), rec->u, rec->v, rec->p);
       return true;
   }

   return false;
}





float material_scattering_pdf(const Inputs* inputs, const int mat_ptr, const ray* r_in, const hit_record* rec, ray* scattered_ray);
float material_scattering_pdf(const Inputs* inputs, const int mat_ptr, const ray* r_in, const hit_record* rec, ray* scattered_ray) {
   float16 mat = inputs->materials[mat_ptr];

   //Lambertian
   if (mat[MATK_TYPE]==0){
       float cosine = dot(rec->normal, normalize(scattered_ray->dir));
       return cosine < 0 ? 0 : cosine / PI;
   }


    return 0;
}























float3 ray_color(const ray* r, ray* bounced_ray, const Inputs* inputs);
float3 ray_color(const ray* r, ray* bounced_ray, const Inputs* inputs) {
    //Test hit on world
    hit_record rec;

    //Hit the background
    if (!worldHit(inputs, inputs->objects, inputs->allocated_objects, r, &rec)){
        return (float3){0,0,0};
    }

    scatter_record srec;
    float3 emitted = material_emitted(inputs, rec.mat_ptr, r, &rec, rec.u, rec.v, rec.p);
    if (!material_scatter(inputs, rec.mat_ptr, r, &rec, &srec)){
        return emitted;
    }

    if (srec.is_specular){
        bounced_ray->orig = srec.specular_ray.orig;
        bounced_ray->dir  = srec.specular_ray.dir;
        return srec.attenuation;
    }

    //Add object pdf to extisting pdf and set it mixed
    int rand_light = rand(inputs->seed) % inputs->allocated_lights;
    pdf_init_object(&srec.scatter_pdf, rand_light, rec.p, inputs);
    srec.scatter_pdf.type = PDFT_MIXTURE;

    //Generate scatter ray
    ray scattered = {rec.p, pdf_generate(&srec.scatter_pdf)};

    //Generate pdf value
    float weight = 1.0f/inputs->allocated_lights;
    float pdf_val = 0.0f;
    for(unsigned int i=0; i<inputs->allocated_lights; i++){
        pdf_val += weight * pdf_value(&srec.scatter_pdf, scattered.dir);
    }

    //Fill bounced ray data
    bounced_ray->orig = scattered.orig;
    bounced_ray->dir  = scattered.dir;
    return emitted
        + srec.attenuation * material_scattering_pdf(inputs, rec.mat_ptr, r, &rec, &scattered) / pdf_val;
}






/**
  *
  */
__kernel void render(uint8 chunk_data,
                     int4 parameters_data,

                     uint8 world_data,
                     __global float16* objects,
                     __global float16* wrapped_objects,
                     __global float16* lights,
                     __global float16* materials,
                     __global float16* textures,

                     __global float16* rays,
                     __global float3* rus,
                     __global int* seeds,

                     __global float3* perlin_ranvec,
                     __global int3* perlin_perms,

                     uint2 iteration_data,
                     __global float4* output) {

   //Get Global 1D ID
   const uint GLOBAL_ID = get_global_id(1)*get_global_size(0)+get_global_id(0);

   //Unpack world data
   const int OBJS_COUNT         = world_data[0];
   const int WRAPPED_OBJS_COUNT = world_data[1];
   const int LIGHTS_COUNT       = world_data[2];
   const int MATS_COUNT         = world_data[3];
   const int TEXS_COUNT         = world_data[4];
   const int ALLOCATED_OBJECTS  = world_data[5];
   const int ALLOCATED_LIGHTS   = world_data[6];

   //Unpack parameters data
   const int SPP                = parameters_data[0];
   const int RANDOM_COUNT       = parameters_data[1];
   const int RUS_COUNT          = parameters_data[2];
   const int PERLIN_POINT_COUNT = parameters_data[3];

   //Setup random generation
   int seed = seeds[GLOBAL_ID % RANDOM_COUNT];

   //Unpack chunk data
   //const uint CHUNK_X     = chunk_data[0];
   //const uint CHUNK_Y     = chunk_data[1];
   const uint CHUNK_W     = chunk_data[2];
   //const uint CHUNK_H     = chunk_data[3];
   //const uint CHUNK_NUM_W = chunk_data[4];
   //const uint CHUNK_NUM_H = chunk_data[5];
   //const uint IMAGE_WIDTH  = CHUNK_W*CHUNK_NUM_W;
   //const uint IMAGE_HEIGHT = CHUNK_H*CHUNK_NUM_H;

   //Unpack iteration data
   const uint ITERATION = iteration_data[0];
   const uint ITERATIONS_COUNT = iteration_data[1];

   //Get the pixel coordinate relative to chunk
   const uint LOCAL_X = get_local_id(0);
   const uint LOCAL_Y = get_local_id(1);
   const uint LOCAL_W = get_local_size(0);
   const uint LOCAL_H = get_local_size(1);
   const uint GROUP_X = get_group_id(0);
   const uint GROUP_Y = get_group_id(1);
   const uint CHUNK_PX = ((GROUP_X * LOCAL_W) + LOCAL_X); //relative to chunk
   const uint CHUNK_PY = ((GROUP_Y * LOCAL_H) + LOCAL_Y); //relative to chunk


   //Pack the input struct
   Inputs inputs;
   inputs.objects = objects;
   inputs.wrapped_objects = wrapped_objects;
   inputs.lights = lights;
   inputs.materials = materials;
   inputs.textures = textures;
   inputs.objs_count = OBJS_COUNT;
   inputs.wrapped_objs_count = WRAPPED_OBJS_COUNT;
   inputs.lights_count = LIGHTS_COUNT;
   inputs.mats_count = MATS_COUNT;
   inputs.texs_count = TEXS_COUNT;

   inputs.allocated_objects = ALLOCATED_OBJECTS;
   inputs.allocated_lights = ALLOCATED_LIGHTS;

   inputs.rus = rus;
   inputs.rus_count = RUS_COUNT;
   inputs.seed = &seed;
   inputs.perlin_ranvec = perlin_ranvec;
   inputs.perlin_perms = perlin_perms;
   inputs.perlin_point_count = PERLIN_POINT_COUNT;

   inputs.GLOBAL_ID = GLOBAL_ID;
   inputs.ITERATION = ITERATION;

   //Final color
   const int i = (CHUNK_PY * CHUNK_W + CHUNK_PX)*SPP;
   for (int s=0;s<SPP;++s){
       inputs.SAMPLE = s;

      //Skip dead rays
      if (rays[i+s][6] != 0){continue;}

      //Unpack the ray from the ray_pool
      ray r;
      r.orig[0] = rays[i+s][0];
      r.orig[1] = rays[i+s][1];
      r.orig[2] = rays[i+s][2];
      r.dir[0]  = rays[i+s][3];
      r.dir[1]  = rays[i+s][4];
      r.dir[2]  = rays[i+s][5];

      //Get the ray_color and retrieve an eventual bounced ray
      ray bounced_ray = {(float3){999,999,999}, (float3){999,999,999}};
      float3 sample_color;
      if (ITERATION == ITERATIONS_COUNT-1){sample_color = (float3){0,0,0};
      }else{sample_color = ray_color(&r, &bounced_ray, &inputs);}

      //Multiply the last sample with the new one. NB: these needs to be initiated to 1.0
      rays[i+s][10] = rays[i+s][10] * sample_color[0];
      rays[i+s][11] = rays[i+s][11] * sample_color[1];
      rays[i+s][12] = rays[i+s][12] * sample_color[2];


      //If I don't need to trace any more ray
      if (bounced_ray.orig[0]==999){
         //Set dead ray
         rays[i+s][6] = 1.0;

         //Save the final ray color into the samples to get averaged
         const int CHUNK_IND = (CHUNK_PY * CHUNK_W) + CHUNK_PX;
         output[CHUNK_IND] += (float4){rays[i+s][10], rays[i+s][11], rays[i+s][12], 1};

         //If the ray is still active
      }else{
         //Unpack the bounced ray into the ray pool for the next pass
         rays[i+s][0] = bounced_ray.orig[0];
         rays[i+s][1] = bounced_ray.orig[1];
         rays[i+s][2] = bounced_ray.orig[2];
         rays[i+s][3] = bounced_ray.dir[0];
         rays[i+s][4] = bounced_ray.dir[1];
         rays[i+s][5] = bounced_ray.dir[2];
      }
   }


   //Save the seed for the next time this kernel gets enqueued.
   seeds[GLOBAL_ID % RANDOM_COUNT] = seed;
}

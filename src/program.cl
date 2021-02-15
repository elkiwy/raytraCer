

/**
 * Rays
 **/

typedef struct ray {
   float3 orig;
   float3 dir;
} ray;

float3 ray_at(const ray* r, float t) {
   return r->orig + t * r->dir;
}




/**
 * Camera
 **/

//typedef struct camera {
//   float3 origin;
//   float3 lower_left_corner;
//   float3 horizontal;
//   float3 vertical;
//
//} camera;
//
//
//void init_camera(camera* c) {
//   float aspect_ratio = 1.0 / 1.0;
//   float viewport_height = 2.0;
//   float viewport_width = aspect_ratio * viewport_height;
//   float focal_length = 1.0;
//
//   c->origin = (float3){0, 0, 0};
//   c->horizontal = (float3){viewport_width, 0.0, 0.0};
//   c->vertical = (float3){0.0, viewport_height, 0.0};
//   c->lower_left_corner = c->origin - c->horizontal/2 - c->vertical/2 - (float3){0, 0, focal_length};
//}
//
//
//ray camera_get_ray(const camera* c, float u, float v) {
//   return (ray){c->origin, c->lower_left_corner + u*c->horizontal + v*c->vertical - c->origin};
//}






/**
 * Hittable
 **/

typedef struct hit_record {
   float3 p;
   float3 normal;
   float t;
   bool front_face;
} hit_record;

void set_face_normal(hit_record* rec, const ray* r, float3 outward_normal) {
   rec->front_face = dot(r->dir, outward_normal) < 0;
   rec->normal = rec->front_face ? outward_normal :-outward_normal;
}

typedef struct hittable {
   float3 bounds1;
   float3 bounds2;
   float3 center;
   float radius;
   uint type;
} hittable;


float length_squared(float3 v){return v.x*v.x + v.y*v.y + v.z*v.z;}
bool hittable_hit(const hittable* obj, const ray* r, double t_min, double t_max, hit_record* rec) {
   if (obj->type==0){
      float3 oc = r->orig - obj->center;
      float a = length_squared(r->dir);
      float half_b = dot(oc, r->dir);
      float c = length_squared(oc) - obj->radius*obj->radius;
      float discriminant = half_b*half_b - a*c;
      if (discriminant < 0) {return false;}

      float sqrtd = sqrt(discriminant);
      float root  = (-half_b - sqrtd) / a;
      if (root < t_min || t_max < root) {
         root = (-half_b + sqrtd) / a;
         if (root < t_min || t_max < root){
            return false;
         }
      }

      rec->t = root;
      rec->p = ray_at(r, rec->t);
      rec->normal = (rec->p - obj->center) / obj->radius;
      float3 outward_normal = (rec->p - obj->center) / obj->radius;
      set_face_normal(rec, r, outward_normal);

      return true;
   }else{

      //TODO
      return false;
   }
}



bool worldHit(hittable* objs, int objs_count, const ray* r, hit_record *rec) {
   hit_record temp_rec;
   bool hitAnything = false;
   float closestSoFar = FLT_MAX;
   for (int i = 0; i < objs_count; i++) {
      if (hittable_hit(&objs[i], r, 0.001, closestSoFar, &temp_rec)) {
         hitAnything = true;
         closestSoFar = temp_rec.t;
         *rec = temp_rec;
      }
   }
   return hitAnything;
}











///Generate random number between 0 and 2147483647
int rand(int* seed) {// 1 <= *seed < m
    int const a = 16807; //ie 7**5
    int const m = 2147483647; //ie 2**31-1
    *seed = (long)(*(seed) * a) % m;
    return(*seed);
}



float3 ray_color(const ray* r, const hittable* objs, int objs_count, __global float3* rus, ray* bounced_ray, int* seed) {
   //Test hit on world
   hit_record rec;
   if (worldHit(objs, objs_count, r, &rec)){
      //Fill bounced ray data
      float3 target = rec.p + rec.normal + rus[rand(seed) % 1024];
      bounced_ray->orig = rec.p;
      bounced_ray->dir  = target - rec.p;

      //Return the color
      return (float3){0,0,0};
   }

   //Background hit
   float3 unit_direction = normalize(r->dir);
   float t = 0.5f * (unit_direction[1] + 1.0);
   return (1.0f-t) * (float3){1.0,1.0,1.0} + t * (float3){0.3, 0.7, 1.0};
}






/**
  *
  */
__kernel void render(uint8 chunk_data,  __global float8* rays, __global float3* rus, __global int* seeds, uint pass_number, __global float4* pixels) {
   //Get Global 1D ID
   const uint GLOBAL_ID = get_global_id(1)*get_global_size(0)+get_global_id(0);

   //TODO: pass SPP from CPU
   const int SPP = 128;

   //Unpack chunk data
   const uint CHUNK_X     = chunk_data[0];
   const uint CHUNK_Y     = chunk_data[1];
   const uint CHUNK_W     = chunk_data[2];
   const uint CHUNK_H     = chunk_data[3];
   const uint CHUNK_NUM_W = chunk_data[4];
   const uint CHUNK_NUM_H = chunk_data[5];
   const uint IMAGE_WIDTH  = CHUNK_W*CHUNK_NUM_W;
   const uint IMAGE_HEIGHT = CHUNK_H*CHUNK_NUM_H;

   //Get the pixel coordinate relative to chunk
   const uint LOCAL_X = get_local_id(0);
   const uint LOCAL_Y = get_local_id(1);
   const uint LOCAL_W = get_local_size(0);
   const uint LOCAL_H = get_local_size(1);
   const uint GROUP_X = get_group_id(0);
   const uint GROUP_Y = get_group_id(1);
   const uint CHUNK_PX = ((GROUP_X * LOCAL_W) + LOCAL_X); //relative to chunk
   const uint CHUNK_PY = ((GROUP_Y * LOCAL_H) + LOCAL_Y); //relative to chunk

   //Setup world
   const int max_objects = 16;
   hittable world[16];
   for(int i=0;i<max_objects;i++){
      world[i].center = (float3){100,100,100};
      world[i].radius = 0.0;
      world[i].type = 0.0;
   }

   world[0].center = (float3){0,0,-1};
   world[0].radius = 0.5;
   world[0].type = 0;
   world[1].center = (float3){0,-100.5,-1};
   world[1].radius = 100;
   world[1].type = 0;

   //Setup random generation
   const int RANDOM_COUNT = CHUNK_W*CHUNK_H*0.05;
   int seed = seeds[GLOBAL_ID % RANDOM_COUNT];

   //Final color
   const int i = (CHUNK_PY * CHUNK_W + CHUNK_PX)*SPP;
   float3 pixel_color = {0,0,0};
   int samples = 0;
   for (int s=0;s<SPP;++s){
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
      pixel_color += ray_color(&r, &world, max_objects, rus, &bounced_ray, &seed);

      //Increment the samples count took for this pixel (NB: SPP - dead rays)
      samples++;

      //Unpack the bounced ray into the ray pool for the next pass
      rays[i+s][0] = bounced_ray.orig[0];
      rays[i+s][1] = bounced_ray.orig[1];
      rays[i+s][2] = bounced_ray.orig[2];
      rays[i+s][3] = bounced_ray.dir[0];
      rays[i+s][4] = bounced_ray.dir[1];
      rays[i+s][5] = bounced_ray.dir[2];
      rays[i+s][6] = bounced_ray.orig[0] == 999; //Set dead ray
   }


   //Save the seed for the next time this kernel gets enqueued.
   seeds[GLOBAL_ID % RANDOM_COUNT] = seed;

   //Write the pixel color to the right absolute pixel of our image
   //(NB: this accumulates the pixel color and samples count to later average those in the CPU)
   const int ABSOLUTE_X = CHUNK_PX + CHUNK_X*CHUNK_W;
   const int ABSOLUTE_Y = CHUNK_PY + CHUNK_Y*CHUNK_H;
   pixels[(ABSOLUTE_Y * IMAGE_WIDTH) + ABSOLUTE_X] += (float4){pixel_color, samples};
}

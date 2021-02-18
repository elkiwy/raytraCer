/**
 * Utils
 **/



//Random in unit sphere to random in hemisphere
float3 rus_to_rhemi(float3 rus, float3 normal){
   if(dot(rus, normal)>0.0){
      return rus;
   }else{
      return -rus;
   }
}

//Random in unit spehere to random in unit vector
float3 rus_to_ruv(float3 rus){
   return normalize(rus);
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
 * Hittable
 **/

struct material;

typedef struct hit_record {
   float3 p;
   float3 normal;
   float t;
   bool front_face;
   struct material* mat_ptr;
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
   struct material* mat;
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
      rec->mat_ptr = obj->mat;

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








/**
 * Materials
 **/

typedef struct material {
   float3 albedo;
   float fuzzy;
   uint type;
   __global float3* rus;
   int rus_count;
   int* seed;
} material;



bool near_zero(float3 v){
   const float s = 1e-8;
   return (fabs(v[0])<s) && (fabs(v[1])<s) && (fabs(v[2])<s);
}

float3 reflect(float3 v, float3 n){
   return v - 2*dot(v, n)*n;
}

bool material_scatter(const material* mat, const ray* r_in, const hit_record* rec, float3* attenuation, ray* scattered_ray) {
   //Lambertian
   if (mat->type==0){
      float3 scatter_direction = rec->normal + rus_to_ruv(mat->rus[rand(mat->seed) % mat->rus_count]);

      if (near_zero(scatter_direction)){
         scatter_direction = rec->normal;
      }

      *scattered_ray = (ray){rec->p, scatter_direction};
      *attenuation = mat->albedo;

      return true;
   //Metal
   }else{

      float3 reflected = reflect(normalize(r_in->dir), rec->normal);
      *scattered_ray = (ray){rec->p, reflected};
      *attenuation = mat->albedo;

      return (dot(scattered_ray->dir, rec->normal) > 0);
   }
}






















float3 ray_color(const ray* r, const hittable* objs, int objs_count, __global float3* rus, ray* bounced_ray, int* seed) {
   //Test hit on world
   hit_record rec;
   if (worldHit(objs, objs_count, r, &rec)){
      ray scattered;
      float3 attenuation;

      if (material_scatter(rec.mat_ptr, r, &rec, &attenuation, &scattered)){
         //Fill bounced ray data
         bounced_ray->orig = scattered.orig;
         bounced_ray->dir  = scattered.dir;
         //Return the attenuation color
         return attenuation;
      }

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
__kernel void render(uint8 chunk_data,
                     int4 parameters_data,
                     __global float16* rays,
                     __global float3* rus,
                     __global int* seeds,
                     uint2 iteration_data,
                     __global float4* output) {

   //Get Global 1D ID
   const uint GLOBAL_ID = get_global_id(1)*get_global_size(0)+get_global_id(0);
   //if (GLOBAL_ID==0){printf("Last pass : %d!\n", pass_number);}


   //Unpack parameters data
   const int SPP          = parameters_data[0];
   const int RANDOM_COUNT = parameters_data[1];
   const int RUS_COUNT    = parameters_data[2];

   //Setup random generation
   int seed = seeds[GLOBAL_ID % RANDOM_COUNT];

   //Unpack chunk data
   const uint CHUNK_X     = chunk_data[0];
   const uint CHUNK_Y     = chunk_data[1];
   const uint CHUNK_W     = chunk_data[2];
   const uint CHUNK_H     = chunk_data[3];
   const uint CHUNK_NUM_W = chunk_data[4];
   const uint CHUNK_NUM_H = chunk_data[5];
   const uint IMAGE_WIDTH  = CHUNK_W*CHUNK_NUM_W;
   const uint IMAGE_HEIGHT = CHUNK_H*CHUNK_NUM_H;

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

   //Setup materials
   const int MAX_OBJECTS = 16;
   material materials[16];
   for(int i=0;i<MAX_OBJECTS;i++){
      materials[i].albedo = (float3){0,0,0};
      materials[i].fuzzy = 0.0;
      materials[i].type = 0;
      materials[i].rus = rus;
      materials[i].rus_count = RUS_COUNT;
      materials[i].seed = &seed;
   }
   materials[0].albedo = (float3){0.8,0.8,0.0};
   materials[0].type = 0;
   materials[1].albedo = (float3){0.7,0.3,0.3};
   materials[1].type = 0;
   materials[2].albedo = (float3){0.8,0.8,0.8};
   materials[2].type = 1;
   materials[3].albedo = (float3){0.8,0.6,0.2};
   materials[3].type = 1;

   //Setup world objects
   hittable world[16];
   for(int i=0;i<MAX_OBJECTS;i++){
      world[i].center = (float3){100,100,100};
      world[i].radius = 0.0;
      world[i].type = 0;
      world[i].mat = 0;
   }
   world[0].center = (float3){0,-100.5,-1};
   world[0].radius = 100;
   world[0].type = 0;
   world[0].mat = &materials[0];
   world[1].center = (float3){0,0,-1};
   world[1].radius = 0.5;
   world[1].type = 0;
   world[1].mat = &materials[1];
   world[2].center = (float3){-1,0,-1};
   world[2].radius = 0.5;
   world[2].type = 0;
   world[2].mat = &materials[2];
   world[3].center = (float3){1,0,-1};
   world[3].radius = 0.5;
   world[3].type = 0;
   world[3].mat = &materials[3];


   //Final color
   const int i = (CHUNK_PY * CHUNK_W + CHUNK_PX)*SPP;
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
      float3 sample_color;
      if (ITERATION == ITERATIONS_COUNT-1){sample_color = (float3){0,0,0};
      }else{sample_color = ray_color(&r, &world, MAX_OBJECTS, rus, &bounced_ray, &seed);}

      //Multiply the last sample with the new one. NB: these needs to be initiated to 1.0
      rays[i+s][10] = rays[i+s][10] * sample_color[0];
      rays[i+s][11] = rays[i+s][11] * sample_color[1];
      rays[i+s][12] = rays[i+s][12] * sample_color[2];


      //If I don't need to trace any more ray
      if (bounced_ray.orig[0]>900){
         //Set dead ray
         rays[i+s][6] = 1.0;

         //Save the final ray color into the samples to get averaged
         const int ABSOLUTE_X = CHUNK_PX + CHUNK_X*CHUNK_W;
         const int ABSOLUTE_Y = CHUNK_PY + CHUNK_Y*CHUNK_H;
         const int ABSOLUTE_IND = (ABSOLUTE_Y * IMAGE_WIDTH) + ABSOLUTE_X;
         output[ABSOLUTE_IND] += (float4){rays[i+s][10], rays[i+s][11], rays[i+s][12], 1};

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

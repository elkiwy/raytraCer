/**
 * Utils
 **/

typedef struct Inputs {
   __global float16* objects;
   __global float16* materials;
   __global float16* textures;
   uint objs_count;
   uint mats_count;
   uint texs_count;

   __global float3* rus;
   int rus_count;
   int* seed;

   __global float3* perlin_ranvec;
   __global int3* perlin_perms;
   int perlin_point_count;
} Inputs;

float3 rus_to_rhemi(float3 rus, float3 normal);
int rand(int* seed);
float rand_float(int* seed);
float3 random_in_unit_sphere(const Inputs* i);
float3 random_in_unit_vector(const Inputs* i);




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

float3 random_in_unit_sphere(const Inputs* i){
   return i->rus[rand(i->seed) % i->rus_count];
}

float3 random_in_unit_vector(const Inputs* i){
   return normalize(random_in_unit_sphere(i));
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





/**
 * Hittable
 **/

//struct material;

typedef struct hit_record {
   float3 p;
   float3 normal;
   float t;
   bool front_face;
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
#define OBJK_POS0 1
#define OBJK_POS1 2
#define OBJK_POS2 3
#define OBJK_RAD  7
#define OBJK_MAT  8
float length_squared(float3 v);
bool hittable_hit(const __global float16* obj_ptr, const ray* r, double t_min, double t_max, hit_record* rec);
bool worldHit(const Inputs* inputs, const ray* r, hit_record *rec);

float length_squared(float3 v){return v.x*v.x + v.y*v.y + v.z*v.z;}




#define PI 3.14159265359
void hittable_get_uv(const __global float16* obj_ptr, float3 p, float* u, float* v);
void hittable_get_uv(const __global float16* obj_ptr, float3 p, float* u, float* v){
   float16 obj = *obj_ptr;
   if(obj[OBJK_TYPE]==0){

      float theta = acos(-p[1]);
      float phi = atan2(-p[2], p[0]) + PI;
      *u = phi / (2*PI);
      *v = theta / PI;

   }else{
      //TODO
   }
}




bool hittable_hit(const __global float16* obj_ptr, const ray* r, double t_min, double t_max, hit_record* rec) {
   float16 obj = *obj_ptr;
   if (obj[OBJK_TYPE]==0){
      float3 center = {obj[OBJK_POS0],obj[OBJK_POS1],obj[OBJK_POS2]};
      float3 oc = r->orig - center;
      float a = length_squared(r->dir);
      float half_b = dot(oc, r->dir);
      float c = length_squared(oc) - obj[OBJK_RAD]*obj[OBJK_RAD];
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
      rec->normal = (rec->p -  center) / obj[OBJK_RAD];
      float3 outward_normal = (rec->p - center) / obj[OBJK_RAD];
      set_face_normal(rec, r, outward_normal);
      hittable_get_uv(obj_ptr, outward_normal, &rec->u, &rec->v);
      rec->mat_ptr = (int)round(obj[OBJK_MAT]);

      return true;
   }else{

      //TODO
      return false;
   }
}





bool worldHit(const Inputs* inputs, const ray* r, hit_record *rec) {
   hit_record temp_rec;
   bool hitAnything = false;
   float closestSoFar = FLT_MAX;
   for (unsigned int i = 0; i < inputs->objs_count; i++) {
      if (hittable_hit(&(inputs->objects[i]), r, 0.001, closestSoFar, &temp_rec)) {
         hitAnything = true;
         closestSoFar = temp_rec.t;
         *rec = temp_rec;
      }
   }
   return hitAnything;
}



















/**
 * Textures
 **/


float noise(const Inputs* inp, float3 p);
float noise(const Inputs* inp, float3 p) {
   float u = p[0] - floor(p[0]);
   float v = p[1] - floor(p[1]);
   float w = p[2] - floor(p[2]);
   //u = u*u*(3-2*u);
   //v = v*v*(3-2*v);
   //w = w*w*(3-2*w);

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

double turb(const Inputs* inp, float3 p, int depth) {
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

      //return (float3){1,1,1};

      //float v = noise(inputs, tex[TEXK_SCALE] * p);
      //return (float3){1,1,1} * 0.5f * (1.0f + v);

      //float v = turb(inputs, tex[TEXK_SCALE] * p, 7);
      //return (float3){1,1,1} * v;

      float v = 0.5 * (1.0f + sin(turb(inputs, p, 7)*10 + tex[TEXK_SCALE]*p[2]));
      return (float3){1,1,1} * v;
   }
   return (float3){0,0,0};
}
















/**
 * Materials
 **/

//typedef struct material {
//   //Material definition
//   uint type;
//   float3 albedo;
//   float fuzz;
//   float ir;
//} material;


#define MATK_TYPE 0
#define MATK_ALB0 1
#define MATK_ALB1 2
#define MATK_ALB2 3
#define MATK_FUZZ 4
#define MATK_IR 5
#define MATK_TEX_PTR 6

bool near_zero(float3 v);
float3 reflect(float3 v, float3 n);
float3 refract(float3 uv, float3 n, float etai_over_etat);
float reflectance(float cosine, float ref_idx) ;
bool material_scatter(const Inputs* inputs, const int mat_ptr, const ray* r_in, const hit_record* rec, float3* attenuation, ray* scattered_ray) ;


bool near_zero(float3 v){
   const float s = 1e-8;
   return (fabs(v[0])<s) && (fabs(v[1])<s) && (fabs(v[2])<s);
}

float3 reflect(float3 v, float3 n){
   return v - 2*dot(v, n)*n;
}

float3 refract(float3 uv, float3 n, float etai_over_etat){
   float cos_theta = fmin((float)dot(-uv, n), 1.0f);
   float3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
   float3 r_out_parallel = -sqrt(fabs(1.0f - length_squared(r_out_perp))) * n;
   return r_out_perp + r_out_parallel;
}


float reflectance(float cosine, float ref_idx) {
   // Use Schlick's approximation for reflectance.
   float r0 = (1.0f-ref_idx) / (1.0f+ref_idx);
   r0 = r0*r0;
   return r0 + (1.0f-r0)*pow((1.0f - cosine),5.0f);
}




bool material_scatter(const Inputs* inputs, const int mat_ptr, const ray* r_in, const hit_record* rec, float3* attenuation, ray* scattered_ray) {
   //Lambertian
   float16 mat = inputs->materials[mat_ptr];
   if (mat[MATK_TYPE]==0){
      float3 scatter_direction = rec->normal + random_in_unit_vector(inputs);

      if (near_zero(scatter_direction)){
         scatter_direction = rec->normal;
      }

      *scattered_ray = (ray){rec->p, scatter_direction};

      //*attenuation = (float3){mat[MATK_ALB0], mat[MATK_ALB1], mat[MATK_ALB2]};
      *attenuation = texture_value(inputs, &(inputs->textures[(int)round(mat[MATK_TEX_PTR])]), rec->u, rec->v, rec->p);

      return true;
   //Metal
   }else if (mat[MATK_TYPE] == 1){

      float3 reflected = reflect(normalize(r_in->dir), rec->normal);
      *scattered_ray = (ray){rec->p, reflected + (mat[MATK_FUZZ] * random_in_unit_sphere(inputs))};
      *attenuation = (float3){mat[MATK_ALB0], mat[MATK_ALB1], mat[MATK_ALB2]};

      return (dot(scattered_ray->dir, rec->normal) > 0);

   //Dielectric
   }else if (mat[MATK_TYPE] == 2){
      *attenuation = (float3){1.0,1.0,1.0};
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

      *scattered_ray = (ray){rec->p, direction};
      return true;
   }

   return false;
}





























float3 ray_color(const ray* r, ray* bounced_ray, const Inputs* inputs);
float3 ray_color(const ray* r, ray* bounced_ray, const Inputs* inputs) {
   //Test hit on world
   hit_record rec;
   if (worldHit(inputs, r, &rec)){
      ray scattered;
      float3 attenuation;

      if (material_scatter(inputs, rec.mat_ptr, r, &rec, &attenuation, &scattered)){
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

                     uint4 world_data,
                     __global float16* objects,
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
   //if (GLOBAL_ID==0){printf("Last pass : %d!\n", pass_number);}




   //Unpack world data
   const int OBJS_COUNT = world_data[0];
   const int MATS_COUNT = world_data[1];
   const int TEXS_COUNT = world_data[2];

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
   inputs.materials = materials;
   inputs.textures = textures;
   inputs.objs_count = OBJS_COUNT;
   inputs.mats_count = MATS_COUNT;
   inputs.texs_count = TEXS_COUNT;
   inputs.rus = rus;
   inputs.rus_count = RUS_COUNT;
   inputs.seed = &seed;
   inputs.perlin_ranvec = perlin_ranvec;
   inputs.perlin_perms = perlin_perms;
   inputs.perlin_point_count = PERLIN_POINT_COUNT;


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
      }else{sample_color = ray_color(&r, &bounced_ray, &inputs);}

      //Multiply the last sample with the new one. NB: these needs to be initiated to 1.0
      rays[i+s][10] = rays[i+s][10] * sample_color[0];
      rays[i+s][11] = rays[i+s][11] * sample_color[1];
      rays[i+s][12] = rays[i+s][12] * sample_color[2];


      //If I don't need to trace any more ray
      if (bounced_ray.orig[0]>900){
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

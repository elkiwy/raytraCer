#define PROGRAM_FILE "src/program.cl"
#define ARRAY_SIZE 1024

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/time.h>

#include <OpenCL/cl.h>

//STB Image write to output png
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "external/stb_image_write.h"



float random_double_unit() {return rand() / (RAND_MAX + 1.0);}
float random_double(float min, float max) {return min + (max-min)*random_double_unit();}
cl_float3 random_in_unit_disk() {
    while (1) {
        cl_float3 p = (cl_float3){{random_double(-1,1), random_double(-1,1), 0}};
        if((p.s[0]*p.s[0] + p.s[1]*p.s[1] + p.s[2]*p.s[2]) < 1){return p;}
    }
}
cl_float3 random_in_unit_sphere(){
    while (1) {
        cl_float3 p = {{random_double(-1,1), random_double(-1,1), random_double(-1,1)}};
        if((p.s[0]*p.s[0] + p.s[1]*p.s[1] + p.s[2]*p.s[2]) < 1){return p;}
    }
}



//Create the OpenCL dedicated device
cl_device_id create_device() {
    cl_platform_id platform;
    cl_device_id dev;

    // Identify a platform
    int err = clGetPlatformIDs(1, &platform, NULL);
    if(err < 0) {perror("Couldn't identify a platform"); exit(1);}

    //Try access the gpu
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &dev, NULL);
    if(err == CL_DEVICE_NOT_FOUND) {
        //Fallback on CPU
        printf("No GPU Found, falling back on CPU.\n");
        err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &dev, NULL);
    }else{printf("GPU Device created succesfully.\n");}
    if(err < 0) {perror("Couldn't access any devices"); exit(1);}
    return dev;
}





//Create program from a file and compile it
cl_program build_program(cl_context ctx, cl_device_id dev, const char* filename) {
    //Read program file and place content into buffer
    FILE* program_handle = fopen(filename, "r");
    if(program_handle == NULL) {perror("Couldn't find the program file"); exit(1);}
    fseek(program_handle, 0, SEEK_END);
    size_t program_size = ftell(program_handle);
    rewind(program_handle);
    char* program_buffer = (char*)malloc(program_size + 1);
    program_buffer[program_size] = '\0';
    fread(program_buffer, sizeof(char), program_size, program_handle);
    fclose(program_handle);

    //Create program from file
    int err;
    cl_program program = clCreateProgramWithSource(ctx, 1, (const char**)&program_buffer, &program_size, &err);
    if(err < 0) {perror("Couldn't create the program"); exit(1);}
    free(program_buffer);

    //Build program
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if(err < 0) {
        // Find size of log and print to std output
        size_t log_size;
        clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
        char* program_log = (char*) malloc(log_size + 1);
        program_log[log_size] = '\0';
        clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG, log_size + 1, program_log, NULL);
        printf("%s\n", program_log);
        free(program_log);
        exit(1);
    }

    return program;
}




float float3_length(cl_float3 v){
    return sqrt((v.s[0] * v.s[0]) + (v.s[1] * v.s[1]) + (v.s[2] * v.s[2]));
}

cl_float3 float3_normalize(cl_float3 v){
    float l = float3_length(v);
    return (cl_float3){{v.s[0]/l, v.s[1]/l, v.s[2]/l}};
}

cl_float3 float3_cross(cl_float3 v1, cl_float3 v2){
    return (cl_float3){{
        v1.s[1] * v2.s[2] - v1.s[2] * v2.s[1],
        v1.s[2] * v2.s[0] - v1.s[0] * v2.s[2],
        v1.s[0] * v2.s[1] - v1.s[1] * v2.s[0]
    }};
}


typedef struct camera {
    float origin[3];
    float lower_left_corner[3];
    float horizontal[3];
    float vertical[3];

    cl_float3 w, u, v;
    float lens_radius;
} camera;


#define PI 3.1415926535897932385
void init_camera(camera* c, cl_float3 from, cl_float3 at, cl_float3 vup, float vfov, float aspect_ratio, float aperture, float focus_dist) {
    float theta = vfov * PI / 180.0;
    float h = tan(theta/2);
    float viewport_height = 2.0*h;
    float viewport_width = aspect_ratio * viewport_height;

    c->w = float3_normalize((cl_float3){{from.s[0]-at.s[0], from.s[1]-at.s[1], from.s[2]-at.s[2]}});
    c->u = float3_normalize(float3_cross(vup, c->w));
    c->v = float3_cross(c->w, c->u);

    float focal_length = 1.0;
    c->origin[0] = from.s[0];
    c->origin[1] = from.s[1];
    c->origin[2] = from.s[2];
    c->horizontal[0] = focus_dist * viewport_width * c->u.s[0];
    c->horizontal[1] = focus_dist * viewport_width * c->u.s[1];
    c->horizontal[2] = focus_dist * viewport_width * c->u.s[2];
    c->vertical[0] = focus_dist * viewport_height * c->v.s[0];
    c->vertical[1] = focus_dist * viewport_height * c->v.s[1];
    c->vertical[2] = focus_dist * viewport_height * c->v.s[2];
    c->lower_left_corner[0] = c->origin[0] - c->horizontal[0]/2.0 - c->vertical[0]/2.0 - focus_dist*c->w.s[0];
    c->lower_left_corner[1] = c->origin[1] - c->horizontal[1]/2.0 - c->vertical[1]/2.0 - focus_dist*c->w.s[1];
    c->lower_left_corner[2] = c->origin[2] - c->horizontal[2]/2.0 - c->vertical[2]/2.0 - focus_dist*c->w.s[2];

    c->lens_radius = aperture / 2.0;
}


cl_float16 get_ray(camera* cam, double s, double t) {
    float offset[3];
    cl_float3 rd = random_in_unit_disk();
    rd.s[0] *= cam->lens_radius;
    rd.s[1] *= cam->lens_radius;
    rd.s[2] *= cam->lens_radius;

    offset[0] = cam->u.s[0]*rd.s[0] + cam->v.s[0]*rd.s[1];
    offset[1] = cam->u.s[1]*rd.s[0] + cam->v.s[1]*rd.s[1];
    offset[2] = cam->u.s[2]*rd.s[0] + cam->v.s[2]*rd.s[1];

    return (cl_float16){{
        /* 0, 1, 2 */
        cam->origin[0] + offset[0],
        cam->origin[1] + offset[1],
        cam->origin[2] + offset[2],
        /* 3, 4, 5 */
        cam->lower_left_corner[0] + s*cam->horizontal[0] + t*cam->vertical[0] - cam->origin[0] - offset[0],
        cam->lower_left_corner[1] + s*cam->horizontal[1] + t*cam->vertical[1] - cam->origin[1] - offset[1],
        cam->lower_left_corner[2] + s*cam->horizontal[2] + t*cam->vertical[2] - cam->origin[2] - offset[2],
        /* 6, 7, 8, 9 */
        0, //dead ray
        0, 0, 0,
        /* 10, 11, 12 */ // RGB
        1, 1, 1,
        /* 13, 14, 15 */
        0, 0, 0
    }};
}
















#define PERLIN_POINT_COUNT 256
void perlin_init(cl_int3* perm_xyz, cl_float3* ranvec){
    for (int i=0;i<PERLIN_POINT_COUNT; ++i){
        ranvec[i] = float3_normalize((cl_float3){{random_double(-1, 1), random_double(-1, 1), random_double(-1, 1)}});
    }
    for(int i=0;i<PERLIN_POINT_COUNT;++i){
        perm_xyz[i].s[0] = i;
        perm_xyz[i].s[1] = i;
        perm_xyz[i].s[2] = i;
    }
    for(int i=0;i<PERLIN_POINT_COUNT;++i){
        for(int j=0;j<3;j++){
            int target = rand() % PERLIN_POINT_COUNT;
            int tmp = perm_xyz[i].s[j];
            perm_xyz[i].s[j] = perm_xyz[target].s[j];
            perm_xyz[target].s[j] = tmp;
        }
    }
}





























typedef enum{
    MAT_LAMBERTIAN = 0,
    MAT_METAL = 1,
    MAT_DIELECTRIC = 2
} material_type;

typedef enum{
    OBJ_SPHERE = 0
} object_type;

typedef enum{
    TEX_SOLID = 0,
    TEX_PERLIN = 1
} texture_type;


int setup_world(cl_float16* objs, cl_float16* mats, cl_float16* texs, unsigned int obj_count, unsigned int mat_count, unsigned int tex_count){
    //Make sure we use clean data
    for(unsigned int i=0;i<obj_count;i++){objs[i] = (cl_float16){{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};}
    for(unsigned int i=0;i<mat_count;i++){mats[i] = (cl_float16){{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};}
    for(unsigned int i=0;i<tex_count;i++){texs[i] = (cl_float16){{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};}

    /* Materials:
    **  0 -> type_flag
    **  1,2,3 -> albedo
    **  4 -> fuzz
    **  5 -> ir
    **  6 -> texture_index
    **/

    /* Objects:
    **  0 -> type_flag
    **  1,2,3 -> point1
    **  4,5,6 -> point2
    **  7 -> radius
    **  8 -> material_index
    **/

    /* Objects:
    **  0 -> type_flag
    **  1,2,3 -> color
    **  4 -> scale
    **/

    int obj_ind = 0;
    int mat_ind = 0;
    int tex_ind = 0;

    texs[tex_ind].s[0] = TEX_SOLID;
    texs[tex_ind].s[1] = 0.25;
    texs[tex_ind].s[2] = 0.25;
    texs[tex_ind].s[3] = 0.25;

    mats[mat_ind].s[0] = MAT_LAMBERTIAN;
    mats[mat_ind].s[6] = tex_ind;

    objs[obj_ind].s[1] = 0;
    objs[obj_ind].s[2] = -1000;
    objs[obj_ind].s[3] = 0;
    objs[obj_ind].s[7] = 1000;
    objs[obj_ind].s[8] = mat_ind;
    mat_ind++;
    obj_ind++;
    tex_ind++;

    /*/
    //Create materials
    for(int a=-5; a<5; a++){
        for(int b=-5; b<5; b++){
            //Type
            objs[obj_ind].s[0] = OBJ_SPHERE;

            //Center
            objs[obj_ind].s[1] = a + 0.9f * random_double_unit();
            objs[obj_ind].s[2] = 0.2f;
            objs[obj_ind].s[3] = b + 0.9f * random_double_unit();

            //Radius
            objs[obj_ind].s[7] = 0.2f;

            double rand_mat = random_double_unit();
            if (rand_mat < 0.8){
                mats[mat_ind].s[0] = MAT_LAMBERTIAN;
                mats[mat_ind].s[1] = random_double_unit() * random_double_unit();
                mats[mat_ind].s[2] = random_double_unit() * random_double_unit();
                mats[mat_ind].s[3] = random_double_unit() * random_double_unit();
                objs[obj_ind].s[8] = mat_ind;
            }else if (rand_mat < 0.95){
                mats[mat_ind].s[0] = MAT_METAL;
                mats[mat_ind].s[1] = (random_double_unit()*0.5f)+0.5f;
                mats[mat_ind].s[2] = (random_double_unit()*0.5f)+0.5f;
                mats[mat_ind].s[3] = (random_double_unit()*0.5f)+0.5f;
                mats[mat_ind].s[4] = random_double_unit()*0.5f;
                objs[obj_ind].s[8] = mat_ind;
            }else{
                mats[mat_ind].s[0] = MAT_DIELECTRIC;
                mats[mat_ind].s[5] = 1.5;
                objs[obj_ind].s[8] = mat_ind;
            }
            mat_ind++;
            obj_ind++;
        }
    }



    mats[mat_ind].s[5] = 1.5;
    mats[mat_ind].s[0] = MAT_DIELECTRIC;
    objs[obj_ind].s[1] = 0;
    objs[obj_ind].s[2] = 1;
    objs[obj_ind].s[3] = 0;
    objs[obj_ind].s[7] = 1.0;
    objs[obj_ind].s[8] = mat_ind;
    mat_ind++;
    obj_ind++;


    mats[mat_ind].s[1] = 0.4;
    mats[mat_ind].s[2] = 0.2;
    mats[mat_ind].s[3] = 0.1;
    mats[mat_ind].s[0] = MAT_LAMBERTIAN;
    objs[obj_ind].s[1] = -4;
    objs[obj_ind].s[2] = 1;
    objs[obj_ind].s[3] = 0;
    objs[obj_ind].s[7] = 1.0;
    objs[obj_ind].s[8] = mat_ind;
    mat_ind++;
    obj_ind++;

    mats[mat_ind].s[1] = 0.7;
    mats[mat_ind].s[2] = 0.6;
    mats[mat_ind].s[3] = 0.5;
    mats[mat_ind].s[4] = 0.1;
    mats[mat_ind].s[0] = MAT_METAL;
    objs[obj_ind].s[1] = 3;
    objs[obj_ind].s[2] = 1;
    objs[obj_ind].s[3] = 0;
    objs[obj_ind].s[7] = 1.0;
    objs[obj_ind].s[8] = mat_ind;
    mat_ind++;
    obj_ind++;

    /*/

    texs[tex_ind].s[0] = TEX_PERLIN;
    texs[tex_ind].s[4] = 10.0;

    mats[mat_ind].s[0] = MAT_LAMBERTIAN;
    mats[mat_ind].s[6] = tex_ind;

    objs[obj_ind].s[1] = 3;
    objs[obj_ind].s[2] = 1;
    objs[obj_ind].s[3] = 0;
    objs[obj_ind].s[7] = 1.0;
    objs[obj_ind].s[8] = mat_ind;
    mat_ind++;
    obj_ind++;
    tex_ind++;
    /**/

    return 1;
}






const unsigned long IMAGE_WIDTH  = 256*2;
const unsigned long IMAGE_HEIGHT = 256*2;

const int CHUNKS_SQRT = 1;
const unsigned long CHUNKS_WIDTH  = IMAGE_WIDTH / CHUNKS_SQRT;
const unsigned long CHUNKS_HEIGHT = IMAGE_HEIGHT / CHUNKS_SQRT;

const int SPP = 128; //Samples per pixels
const int ITERATIONS = 16; //Ray recursion count

int main() {
    /**
    ** OpenCL Configuration
     **/
    //Update random seed
    struct timeval time; gettimeofday(&time, NULL); srand(time.tv_usec);

    // Create device and context, and build program
    cl_int err;
    cl_device_id device = create_device();
    cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
    if(err < 0) {perror("Couldn't create a context"); exit(1);}
    cl_program program = build_program(context, device, PROGRAM_FILE);

    //Setup parallelization parameters
    size_t global_size[2] = {CHUNKS_WIDTH, CHUNKS_HEIGHT};
    size_t local_size[2] = {16, 16}; //MAX_WORK_GROUP_SIZE = 256, so max local size is 16*16

    //Queue and kernel
    cl_command_queue queue = clCreateCommandQueue(context, device, 0, &err);
    if(err < 0) {perror("Couldn't create a command queue"); exit(1);};
    cl_kernel kernel = clCreateKernel(program, "render", &err); // <- kernel function name
    if(err < 0) {perror("Couldn't create a kernel"); exit(1);};


    /**
    ** Kernel arguments setup
     **/
    int argc = 0;
    const int F_R_C  = CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR;
    const int F_RW_C = CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR;
    const int F_RW_U = CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR;

    //Chunk_data = {x, y, w, h, num_h, num_v, ?, ?}
    cl_uint8 chunk_data = {{0,0,CHUNKS_WIDTH,CHUNKS_HEIGHT, CHUNKS_SQRT, CHUNKS_SQRT, 0, 0}};
    printf("IMAGE_WIDTH: %lu\n", IMAGE_WIDTH);
    printf("IMAGE_HEIGHT: %lu\n", IMAGE_HEIGHT);
    printf("CHUNKS_WIDTH: %lu\n", CHUNKS_WIDTH);
    printf("CHUNKS_HEIGHT: %lu\n", CHUNKS_HEIGHT);
    printf("CHUNKS NUMBER: %d\n", CHUNKS_SQRT * CHUNKS_SQRT);
    err = clSetKernelArg(kernel, argc, sizeof(cl_uint8), &chunk_data);
    const int ARGC_CHUNK_DATA = argc;
    argc++;

    //Parameters_data = {SPP, RANDOM_COUNT}
    const int RANDOM_SEEDS_COUNT = CHUNKS_HEIGHT*CHUNKS_WIDTH*0.5;
    const int RUS_COUNT = 1024*4;
    cl_int4 parameters_data = {{SPP, RANDOM_SEEDS_COUNT, RUS_COUNT, PERLIN_POINT_COUNT}};
    err = clSetKernelArg(kernel, argc, sizeof(cl_int4), &parameters_data);
    const int ARGC_PARAMETERS_DATA = argc;
    argc++;



    //World setup
    const int OBJS_COUNT = 128;
    cl_float16 objs[OBJS_COUNT];
    const int MATS_COUNT = 128;
    cl_float16 mats[MATS_COUNT];
    const int TEXS_COUNT = 128;
    cl_float16 texs[TEXS_COUNT];
    err = setup_world(&objs[0], &mats[0], &texs[0], OBJS_COUNT, MATS_COUNT, TEXS_COUNT);
    cl_uint4 world_data = {{OBJS_COUNT, MATS_COUNT, TEXS_COUNT}};
    err = clSetKernelArg(kernel, argc, sizeof(cl_uint4), &world_data);
    argc++;

    //Objects
    cl_mem objects_buffer = clCreateBuffer(context, F_R_C, OBJS_COUNT * sizeof(cl_float16), objs, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &objects_buffer);
    argc++;

    //Materials
    cl_mem materials_buffer = clCreateBuffer(context, F_R_C, MATS_COUNT * sizeof(cl_float16), mats, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &materials_buffer);
    argc++;

    //Textures
    cl_mem textures_buffer = clCreateBuffer(context, F_R_C, TEXS_COUNT * sizeof(cl_float16), texs, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &textures_buffer);
    argc++;




    //Ray Pools = {0,1,2 -> Origin; 3,4,5 -> Direction; 6 -> Used flag; 7 -> ??}
    cl_float3 from = {{13, 2, 3}};
    cl_float3   to = {{ 0, 0, 0}};
    cl_float3  vup = {{ 0, 1, 0}};
    float dist_to_focus = 10.0;
    float aperture = 0.1;
    camera cam; init_camera(&cam, from, to, vup, 20.0, 1.0, aperture, dist_to_focus);
    const unsigned long RAY_POOL_SIZE = (CHUNKS_WIDTH * CHUNKS_HEIGHT) * SPP;
    printf("Creating ray_pool of size %luk\n", RAY_POOL_SIZE/1000);fflush(stdout);
    cl_float16* ray_pool = malloc(sizeof(cl_float16)*RAY_POOL_SIZE);
    cl_mem ray_pool_buffer  = clCreateBuffer(context, F_RW_U, RAY_POOL_SIZE * sizeof(cl_float16), ray_pool, &err);
    if(err < 0) {perror("Couldn't create raypool a buffer"); exit(1);};
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem),   &ray_pool_buffer);
    printf("Created ray_pool\n");fflush(stdout);
    argc++;

    //Random in unit sphere
    cl_float3 rus[RUS_COUNT];
    for (int i=0; i<RUS_COUNT; ++i){rus[i] = random_in_unit_sphere();}
    cl_mem rus_buffer = clCreateBuffer(context, F_R_C, RUS_COUNT * sizeof(cl_float3), rus, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &rus_buffer);
    argc++;

    //Random seeds
    cl_int random_seeds[RANDOM_SEEDS_COUNT];
    for (int i=0; i<RANDOM_SEEDS_COUNT; ++i){random_seeds[i] = rand();}
    cl_mem random_seeds_buffer = clCreateBuffer(context, F_R_C, RANDOM_SEEDS_COUNT * sizeof(cl_int), random_seeds, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &random_seeds_buffer);
    argc++;


    //Perlin
    cl_int3 perm_xyz[PERLIN_POINT_COUNT];
    cl_float3 ranvec[PERLIN_POINT_COUNT];
    perlin_init(&perm_xyz[0], &ranvec[0]);
    cl_mem perlin_ranvec_buffer = clCreateBuffer(context, F_R_C, PERLIN_POINT_COUNT * sizeof(cl_float3), ranvec, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &perlin_ranvec_buffer);
    argc++;
    cl_mem perlin_perm_buffer = clCreateBuffer(context, F_R_C, PERLIN_POINT_COUNT * sizeof(cl_int3), perm_xyz, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &perlin_perm_buffer);
    argc++;



    //Iterations data
    cl_uint2 iteration_data = (cl_uint2){{0, ITERATIONS}};
    err = clSetKernelArg(kernel, argc, sizeof(cl_uint2),  &iteration_data);
    const int ARGC_ITERATIONS_DATA = argc;
    argc++;

    //Pixels
    const unsigned long PIXEL_COUNT_PER_CHUNK = CHUNKS_WIDTH * CHUNKS_HEIGHT;
    printf("Creating output of size %luk\n", PIXEL_COUNT_PER_CHUNK/1000);fflush(stdout);
    cl_float4 output[PIXEL_COUNT_PER_CHUNK];
    cl_mem output_buffer = clCreateBuffer(context, F_RW_C, PIXEL_COUNT_PER_CHUNK * sizeof(cl_float4), output, &err);
    if(err < 0) {perror("Couldn't create pixels a buffer"); exit(1);};
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &output_buffer);
    printf("Created output\n");fflush(stdout);


    printf("Kernel ready for execution:\n");fflush(stdout);




    /**
    ** Pass execution
     **/
    const unsigned long PIXEL_COUNT_PER_IMAGE = IMAGE_WIDTH * IMAGE_HEIGHT;
    const int CHANNELS = 3;
    unsigned char* pixelBytes = malloc(sizeof(unsigned char)*(PIXEL_COUNT_PER_IMAGE*3));


    for (int chunk_y = 0; chunk_y < CHUNKS_SQRT; chunk_y++){
        for (int chunk_x = 0; chunk_x < CHUNKS_SQRT; chunk_x++){
            //Update chunk data
            chunk_data.s[0] = chunk_x;
            chunk_data.s[1] = chunk_y;
            clSetKernelArg(kernel, ARGC_CHUNK_DATA, sizeof(cl_uint8), &chunk_data); // pass

            //Update ray pool data
            const float INV_W = 1.0f/(IMAGE_WIDTH-1);
            const float INV_H = 1.0f/(IMAGE_HEIGHT-1);
            for (unsigned long py=0; py<CHUNKS_HEIGHT; ++py){
                for (unsigned long px=0; px<CHUNKS_WIDTH; ++px){
                    for (int s=0;s<SPP;++s){
                        float u = ((float)(px + CHUNKS_WIDTH*chunk_x) + random_double_unit()) * INV_W;
                        float v = ((float)(py + CHUNKS_HEIGHT*chunk_y) + random_double_unit()) * INV_H;
                        int ray_index = ((py * CHUNKS_WIDTH + px) * SPP) + s;
                        ray_pool[ray_index] = get_ray(&cam, u, v);
                    }
                }
            }

            //Update output buffer
            for (unsigned long i=0;i<PIXEL_COUNT_PER_CHUNK;i++){output[i].s[0] = 0; output[i].s[1] = 0; output[i].s[2] = 0; output[i].s[3] = 0;}

            //Ensure the buffer gets updated with the new set of ray each iteration (not necessary on certain GPU, but needed for others)
            err = clEnqueueWriteBuffer(queue, ray_pool_buffer, CL_FALSE, 0, RAY_POOL_SIZE * sizeof(cl_float16), ray_pool, 0, NULL, NULL);
            err = clEnqueueWriteBuffer(queue, output_buffer, CL_FALSE, 0, PIXEL_COUNT_PER_CHUNK * sizeof(cl_float4), output, 0, NULL, NULL);

            for (int i=0; i<ITERATIONS; ++i){
                //Update pass number
                iteration_data = (cl_uint2){{i, ITERATIONS}};
                clSetKernelArg(kernel, ARGC_ITERATIONS_DATA, sizeof(cl_uint2), &iteration_data);

                // Deploy the kernel to a device.
                err = clEnqueueNDRangeKernel(queue, kernel, 2, NULL, global_size, local_size, 0, NULL, NULL);
                //Synchronize en
                err = clEnqueueReadBuffer(queue, output_buffer, CL_TRUE, 0, sizeof(output), output, 0, NULL, NULL);
                //printf("Chunk %d %d pass %d done.\n", chunk_x, chunk_y, iteration_data.s[0]);
            }
            printf("Chunk %d %d processed, now writing...\n", chunk_x, chunk_y);



            /**
            ** Kernel output read and processing
            **/
            //Process the output into an image
            for (unsigned long i=0; i<PIXEL_COUNT_PER_CHUNK; i++){
                const unsigned long current_chunk_x = i % CHUNKS_WIDTH;
                const unsigned long current_chunk_y = i / CHUNKS_WIDTH;
                const unsigned long image_x = current_chunk_x + (chunk_x*CHUNKS_WIDTH);
                const unsigned long image_y = current_chunk_y + (chunk_y*CHUNKS_HEIGHT);
                const unsigned long image_ind = (image_y * IMAGE_WIDTH) + image_x;

                if (i==0){printf("\n\n%d pixel iteration %d: %f %f %f (%f samples)\n", i, 0, output[i].s[0], output[i].s[1], output[i].s[2], output[i].s[3]);}
                //Average of samples
                float samples_for_iteration = output[i].s[3];
                float r = output[i].s[0]/samples_for_iteration;
                float g = output[i].s[1]/samples_for_iteration;
                float b = output[i].s[2]/samples_for_iteration;

                //Gamma correction
                r = sqrt(r);
                g = sqrt(g);
                b = sqrt(b);
                if (i==0){printf("%d final values : %f %f %f \n", i, r, g, b);}

                //Write to pixel buffer
                pixelBytes[(image_ind*3)+0] = (unsigned char)(r*255);
                pixelBytes[(image_ind*3)+1] = (unsigned char)(g*255);
                pixelBytes[(image_ind*3)+2] = (unsigned char)(b*255);
                //if (i>180220){printf("done! %d\n", image_ind);fflush(stdout);}
            }
            printf("Chunk %d %d completed.\n", chunk_x, chunk_y);


        }
    }




    //Flip the final image and save it
    stbi_flip_vertically_on_write(1);
    stbi_write_png("test_opencl.png", IMAGE_WIDTH, IMAGE_HEIGHT, CHANNELS, pixelBytes, IMAGE_WIDTH * CHANNELS);


    //Deallocate resources
    free(ray_pool);
    free(pixelBytes);
    clReleaseKernel(kernel);

    clReleaseMemObject(perlin_ranvec_buffer);
    clReleaseMemObject(perlin_perm_buffer);

    clReleaseMemObject(objects_buffer);
    clReleaseMemObject(materials_buffer);
    clReleaseMemObject(textures_buffer);
    clReleaseMemObject(output_buffer);
    clReleaseMemObject(ray_pool_buffer);
    clReleaseMemObject(random_seeds_buffer);
    clReleaseMemObject(rus_buffer);
    clReleaseCommandQueue(queue);
    clReleaseProgram(program);
    clReleaseContext(context);
    return 0;
}

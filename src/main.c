#define PROGRAM_FILE "src/program.cl"
#define ARRAY_SIZE 1024

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/time.h>

#include <OpenCL/cl.h>

//Include OpemMP for multithreading
#include <omp.h>


//STB Image write to output png
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "external/stb_image_write.h"



#include "utils.h"
#include "camera.h"
#include "objects.h"




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















































int setup_world(cl_float16* objs, cl_float16* wrapped_objs, cl_float16* lights,  cl_float16* mats, cl_float16* texs,
                unsigned int obj_count, unsigned int wrapped_obj_count, unsigned int lights_count, unsigned int mat_count, unsigned int tex_count,
                int* allocated_objects,
                int* allocated_lights){

    //Make sure we use clean data
    for(unsigned int i=0;i<obj_count;i++){objs[i] = (cl_float16){{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};}
    for(unsigned int i=0;i<mat_count;i++){mats[i] = (cl_float16){{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};}
    for(unsigned int i=0;i<tex_count;i++){texs[i] = (cl_float16){{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};}
    for(unsigned int i=0;i<lights_count;i++){lights[i] = (cl_float16){{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};}
    for(unsigned int i=0;i<wrapped_obj_count;i++){wrapped_objs[i] = (cl_float16){{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};}


    int wrapped_obj_ind = 0;
    int obj_ind = 0;
    int light_ind = 0;
    int mat_ind = 0;
    int tex_ind = 0;

    int red_tex = add_to_list(make_texture_solid(0.65, 0.05, 0.05), texs, &tex_ind);
    int mat_red = add_to_list(make_material_lambertian(red_tex), mats, &mat_ind);

    int white_tex = add_to_list(make_texture_solid(0.73, 0.73, 0.73), texs, &tex_ind);
    int mat_white = add_to_list(make_material_lambertian(white_tex), mats, &mat_ind);

    int green_tex = add_to_list(make_texture_solid(0.12, 0.45, 0.15), texs, &tex_ind);
    int mat_green = add_to_list(make_material_lambertian(green_tex), mats, &mat_ind);

    int light_tex = add_to_list(make_texture_solid(15, 15, 15), texs, &tex_ind);
    int mat_light = add_to_list(make_material_light(light_tex), mats, &mat_ind);

    int isotropic_tex = add_to_list(make_texture_solid(1, 1, 1), texs, &tex_ind);
    int mat_isotropic = add_to_list(make_material_isotropic(isotropic_tex), mats, &mat_ind);

    int mat_glass = add_to_list(make_material_glass(1.5), mats, &mat_ind);
    mat_ind++;


    int metal_tex = add_to_list(make_texture_solid(0.8, 0.85, 0.88), texs, &tex_ind);
    int mat_metal = add_to_list(make_material_metal(0, metal_tex), mats, &mat_ind);




    add_to_list(make_rect(0,0, 555,555, 555, AXIS_YZ, mat_green), objs, &obj_ind);
    add_to_list(make_rect(0,0, 555,555,   0, AXIS_YZ, mat_red), objs, &obj_ind);

    cl_float16 light = make_rect(213,227, 343,332, 554, AXIS_XZ, mat_light);
    add_to_list(light, lights, &light_ind);
    light = make_flip_face(wrapped_objs, &wrapped_obj_ind, light);
    add_to_list(light, objs, &obj_ind);

    add_to_list(make_rect(0,0, 555,555,   0, AXIS_XZ, mat_white), objs, &obj_ind);
    add_to_list(make_rect(0,0, 555,555, 555, AXIS_XZ, mat_white), objs, &obj_ind);
    add_to_list(make_rect(0,0, 555,555, 555, AXIS_XY, mat_white), objs, &obj_ind);



    //cl_float16 box1 = make_box(wrapped_objs, &wrapped_obj_ind, 0,0,0,  165,330,165, mat_white);
    //box1 = make_rotated(wrapped_objs, &wrapped_obj_ind, box1, AXIS_Y, 15);
    //box1 = make_translated(wrapped_objs, &wrapped_obj_ind, box1, 265,0,295);
    //add_object(box1, objs, &obj_ind);

    //cl_float16 box2 = make_box(wrapped_objs, &wrapped_obj_ind, 0,0,0,  165,165,165, mat_white);
    //box2 = make_rotated(wrapped_objs, &wrapped_obj_ind, box2, AXIS_Y, -18);
    //box2 = make_translated(wrapped_objs, &wrapped_obj_ind, box2, 130,0,65);
    //add_object(box2, objs, &obj_ind);

    //add_object(make_constant_medium_sphere(    275, 300, 275, 2500,  0.0005, mat_isotropic), objs, &obj_ind);


    *allocated_objects = obj_ind;
    *allocated_lights = light_ind;
    return 1;
}








int get_optimal_chunk_splitting(const int img_w, const int img_h, const int samples_per_pixel, cl_device_id dev){
    //unsigned long long* value;
    //size_t valueSize;
    //clGetDeviceInfo(dev, CL_DEVICE_MAX_MEM_ALLOC_SIZE, 0, NULL, &valueSize);
    //value = (unsigned long long*) malloc(valueSize);
    //clGetDeviceInfo(dev, CL_DEVICE_MAX_MEM_ALLOC_SIZE, valueSize, value, NULL);
    //printf("Device: %llu\n", *value);
    //free(value);
    //unsigned long DEVICE_MAX_MEM_ALLOC = *value;
    unsigned long DEVICE_MAX_MEM_ALLOC = 512*1024*1024; //Override with 512 MB

    const int ptr_size = 8;
    int tmp_sqrt = 1;
    unsigned long tmp_w = img_w / tmp_sqrt;
    unsigned long tmp_h = img_h / tmp_sqrt;
    unsigned long tmp_ray_count = (tmp_w * tmp_h) * samples_per_pixel;
    unsigned long tmp_pool_bytes = tmp_ray_count*ptr_size;

    while(tmp_pool_bytes >= DEVICE_MAX_MEM_ALLOC){
        tmp_sqrt *= 2;
        printf("%lu (%lu MB) rays are too much for device (MAX_MEM: %lu), splitting in %d chunks...\n", tmp_ray_count, tmp_pool_bytes/(1024*1024), DEVICE_MAX_MEM_ALLOC/(1024*1024), tmp_sqrt);fflush(stdout);
        tmp_w = img_w / tmp_sqrt;
        tmp_h = img_h / tmp_sqrt;
        tmp_ray_count = (tmp_w * tmp_h) * samples_per_pixel;
        tmp_pool_bytes = tmp_ray_count*ptr_size;
    }

    printf("%lu (%lu MB) rays are fine for device (MAX_MEM: %lu), splitted in %d chunks.\n", tmp_ray_count, tmp_pool_bytes/(1024*1024), DEVICE_MAX_MEM_ALLOC/(1024*1024), tmp_sqrt);fflush(stdout);
    return tmp_sqrt;
}















typedef struct scene_settings{
    int res_width;
    int res_height;
    int samples_per_pixel;
    int iterations;
}scene_settings;

typedef struct scene{
    cl_float16* objs;
    cl_float16* lights;
    cl_float16* wrapped_objs;
    cl_float16* mats;
    cl_float16* texs;

    int objs_max;
    int lights_max;
    int wrapped_objs_max;
    int mats_max;
    int texs_max;

    int objs_count;
    int lights_count;
}scene;


void render_scene(const scene* s, const scene_settings* settings){
    cl_int err;

    //Unpack settings
    const unsigned long IMAGE_WIDTH  = settings->res_width;
    const unsigned long IMAGE_HEIGHT = settings->res_height;
    const int SPP = settings->samples_per_pixel;
    const int ITERATIONS = settings->iterations;

    //Unpack scene
    const int OBJS_COUNT = s->objs_max;
    cl_float16* objs = s->objs;
    int allocated_objects = s->objs_count;
    const int WRAPPED_OBJS_COUNT = s->wrapped_objs_max;
    cl_float16* wrapped_objs = s->wrapped_objs;
    const int LIGHTS_COUNT = s->lights_max;
    cl_float16* lights = s->lights;
    int allocated_lights = s->lights_count;
    const int MATS_COUNT = s->mats_max;
    cl_float16* mats = s->mats;
    const int TEXS_COUNT = s->texs_max;
    cl_float16* texs = s->texs;



    /**
    ** OpenCL Configuration
     **/
    //Update random seed
    struct timeval time; gettimeofday(&time, NULL); srand(time.tv_usec);
    printf("seed: %d\n",time.tv_usec);

    // Create device and context, and build program
    cl_device_id device = create_device();
    cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
    if(err < 0) {perror("Couldn't create a context"); exit(1);}
    cl_program program = build_program(context, device, PROGRAM_FILE);


    const int CHUNKS_SQRT = get_optimal_chunk_splitting(IMAGE_WIDTH, IMAGE_HEIGHT, SPP, device);
    const unsigned long CHUNKS_WIDTH  = IMAGE_WIDTH / CHUNKS_SQRT;
    const unsigned long CHUNKS_HEIGHT = IMAGE_HEIGHT / CHUNKS_SQRT;


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
    if(err != 0) {perror("Couldn't set kernel arg chunk_data"); exit(1);};
    const int ARGC_CHUNK_DATA = argc;
    argc++;

    //Parameters_data = {SPP, RANDOM_COUNT}
    const int RANDOM_SEEDS_COUNT = CHUNKS_HEIGHT*CHUNKS_WIDTH*0.5;
    const int RUS_COUNT = 1024*4;
    cl_int4 parameters_data = {{SPP, RANDOM_SEEDS_COUNT, RUS_COUNT, PERLIN_POINT_COUNT}};
    err = clSetKernelArg(kernel, argc, sizeof(cl_int4), &parameters_data);
    if(err != 0) {perror("Couldn't set kernel arg params_data"); exit(1);};
    const int ARGC_PARAMETERS_DATA = argc;
    argc++;



    //World setup
    cl_uint8 world_data = {{OBJS_COUNT, WRAPPED_OBJS_COUNT, LIGHTS_COUNT, MATS_COUNT, TEXS_COUNT, allocated_objects, allocated_lights}};
    err = clSetKernelArg(kernel, argc, sizeof(cl_uint8), &world_data);
    if(err != 0) {perror("Couldn't set kernel arg world_data"); exit(1);};
    argc++;

    //Objects
    cl_mem objects_buffer = clCreateBuffer(context, F_R_C, OBJS_COUNT * sizeof(cl_float16), objs, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &objects_buffer);
    if(err != 0) {perror("Couldn't set kernel arg objects"); exit(1);};
    argc++;

    //wrapped Objects
    cl_mem wrapped_objects_buffer = clCreateBuffer(context, F_R_C, WRAPPED_OBJS_COUNT * sizeof(cl_float16), wrapped_objs, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &wrapped_objects_buffer);
    if(err != 0) {perror("Couldn't set kernel arg wrapped_objs"); exit(1);};
    argc++;

    //Lights
    cl_mem lights_buffer = clCreateBuffer(context, F_R_C, LIGHTS_COUNT * sizeof(cl_float16), lights, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &lights_buffer);
    if(err != 0) {perror("Couldn't set kernel arg lights"); exit(1);};
    argc++;

    //Materials
    cl_mem materials_buffer = clCreateBuffer(context, F_R_C, MATS_COUNT * sizeof(cl_float16), mats, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &materials_buffer);
    if(err != 0) {perror("Couldn't set kernel arg materials"); exit(1);};
    argc++;

    //Textures
    cl_mem textures_buffer = clCreateBuffer(context, F_R_C, TEXS_COUNT * sizeof(cl_float16), texs, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &textures_buffer);
    if(err != 0) {perror("Couldn't set kernel arg textures"); exit(1);};
    argc++;




    // 1024 -> float3
    // 256 -> float3
    // 256 -> int3
    //
    // 256 * 256 -> float4           =     1.000.000 bytes
    // 959                             4.022.337.536
    // 960 * 256 * 256 -> float16    = 4.026.531.840 byte


    //Ray Pools = {0,1,2 -> Origin; 3,4,5 -> Direction; 6 -> Used flag; 7 -> ??}
    cl_float3 from = {{278, 278, -800}};
    cl_float3   to = {{278, 278, 0}};
    cl_float3  vup = {{ 0, 1, 0}};
    float dist_to_focus = 10.0;
    float aperture = 0.0;
    float vfov = 40.0;
    camera cam; init_camera(&cam, from, to, vup, vfov, 1.0, aperture, dist_to_focus);

    const unsigned long RAY_POOL_SIZE = (CHUNKS_WIDTH * CHUNKS_HEIGHT) * SPP;
    printf("Creating ray_pool of size %lu MB bytes\n", (8*RAY_POOL_SIZE/(1024*1024)));fflush(stdout);
    cl_float16* ray_pool = malloc(sizeof(cl_float16)*RAY_POOL_SIZE);
    cl_mem ray_pool_buffer  = clCreateBuffer(context, F_RW_U, RAY_POOL_SIZE * sizeof(cl_float16), ray_pool, &err);
    if(err < 0) {perror("Couldn't create raypool a buffer"); exit(1);};
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem),   &ray_pool_buffer);
    if(err != 0) {perror("Couldn't set kernel arg rays"); exit(1);};
    printf("Created ray_pool\n");fflush(stdout);
    argc++;

    //Random in unit sphere
    cl_float3 rus[RUS_COUNT];
    for (int i=0; i<RUS_COUNT; ++i){rus[i] = random_in_unit_sphere();}
    cl_mem rus_buffer = clCreateBuffer(context, F_R_C, RUS_COUNT * sizeof(cl_float3), rus, &err);
    if(err < 0) {perror("Couldn't create rus buffer"); exit(1);};
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &rus_buffer);
    if(err < 0) {perror("Couldn't set kernel arg for rus buffer"); exit(1);};
    argc++;

    //Random seeds
    cl_int random_seeds[RANDOM_SEEDS_COUNT];
    for (int i=0; i<RANDOM_SEEDS_COUNT; ++i){random_seeds[i] = rand();}
    cl_mem random_seeds_buffer = clCreateBuffer(context, F_R_C, RANDOM_SEEDS_COUNT * sizeof(cl_int), random_seeds, &err);
    if(err < 0) {perror("Couldn't create seeds buffer"); exit(1);};
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &random_seeds_buffer);
    if(err < 0) {perror("Couldn't set kernel arg for seeds buffer"); exit(1);};
    argc++;


    //Perlin
    cl_int3 perm_xyz[PERLIN_POINT_COUNT];
    cl_float3 ranvec[PERLIN_POINT_COUNT];
    perlin_init(&perm_xyz[0], &ranvec[0]);
    cl_mem perlin_ranvec_buffer = clCreateBuffer(context, F_R_C, PERLIN_POINT_COUNT * sizeof(cl_float3), ranvec, &err);
    if(err < 0) {perror("Couldn't create ranvec buffer"); exit(1);};
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &perlin_ranvec_buffer);
    if(err < 0) {perror("Couldn't set kernel arg for ranvec buffer"); exit(1);};
    argc++;
    cl_mem perlin_perm_buffer = clCreateBuffer(context, F_R_C, PERLIN_POINT_COUNT * sizeof(cl_int3), perm_xyz, &err);
    if(err < 0) {perror("Couldn't create perms buffer"); exit(1);};
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &perlin_perm_buffer);
    if(err < 0) {perror("Couldn't set kernel arg for perms buffer"); exit(1);};
    argc++;



    //Iterations data
    cl_uint2 iteration_data = (cl_uint2){{0, ITERATIONS}};
    err = clSetKernelArg(kernel, argc, sizeof(cl_uint2),  &iteration_data);
    if(err < 0) {perror("Couldn't set kernel arg for iteration_data buffer"); exit(1);};
    const int ARGC_ITERATIONS_DATA = argc;
    argc++;

    //Pixels
    const unsigned long PIXEL_COUNT_PER_CHUNK = CHUNKS_WIDTH * CHUNKS_HEIGHT;
    printf("Creating output of size %luk\n", PIXEL_COUNT_PER_CHUNK/1000);fflush(stdout);
    cl_float4 output[PIXEL_COUNT_PER_CHUNK];
    cl_mem output_buffer = clCreateBuffer(context, F_RW_C, PIXEL_COUNT_PER_CHUNK * sizeof(cl_float4), output, &err);
    if(err < 0) {perror("Couldn't create pixels a buffer"); exit(1);};
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &output_buffer);
    if(err < 0) {perror("Couldn't set kernel arg for pixels a buffer"); exit(1);};
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

            const unsigned long THREADS = 8;
            const unsigned long step = CHUNKS_HEIGHT/THREADS;

            const unsigned long CHUNK_OFFX = CHUNKS_WIDTH*chunk_x;
            const unsigned long CHUNK_OFFY = CHUNKS_HEIGHT*chunk_y;

            printf("\nGenerating %luk rays...\n", CHUNKS_HEIGHT*CHUNKS_WIDTH*SPP/1000);
            #pragma omp parallel for
            for(unsigned int k=0; k<THREADS; ++k){
                for (int sy=step-1; sy>=0; --sy){
                    unsigned long py = sy + k*step;

                    for (unsigned long px=0; px<CHUNKS_WIDTH; ++px){
                        for (int s=0;s<SPP;++s){
                            float u = ((float)(px + CHUNK_OFFX) + random_double_unit()) * INV_W;
                            float v = ((float)(py + CHUNK_OFFY) + random_double_unit()) * INV_H;
                            int ray_index = ((py * CHUNKS_WIDTH + px) * SPP) + s;
                            get_ray(&cam, u, v, &ray_pool[ray_index]);
                        }
                    }
                }
                //int thread_id = omp_get_thread_num();
                //printf("Thread %d done.\n", thread_id);
            }
            printf("Done\n");


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
                printf("Chunk %d %d pass %d done.\n", chunk_x, chunk_y, iteration_data.s[0]);
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

                //if (i==0){printf("\n\n%d pixel iteration %d: %f %f %f (%f samples)\n", i, 0, output[i].s[0], output[i].s[1], output[i].s[2], output[i].s[3]);}
                //Average of samples
                float inv_samples_for_iteration = 1.0f/output[i].s[3];
                float r = output[i].s[0]*inv_samples_for_iteration;
                float g = output[i].s[1]*inv_samples_for_iteration;
                float b = output[i].s[2]*inv_samples_for_iteration;


                // Replace NaN components with zero. See explanation in Ray Tracing: The Rest of Your Life.
                if (r != r) r = 0.0;
                if (g != g) g = 0.0;
                if (b != b) b = 0.0;

                //Gamma correction
                r = fmax(fmin(sqrt(r), 0.99), 0.0);
                g = fmax(fmin(sqrt(g), 0.99), 0.0);
                b = fmax(fmin(sqrt(b), 0.99), 0.0);

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
    clReleaseMemObject(wrapped_objects_buffer);
    clReleaseMemObject(materials_buffer);
    clReleaseMemObject(textures_buffer);
    clReleaseMemObject(output_buffer);
    clReleaseMemObject(ray_pool_buffer);

    clReleaseMemObject(random_seeds_buffer);
    clReleaseMemObject(rus_buffer);

    clReleaseCommandQueue(queue);
    clReleaseProgram(program);
    clReleaseContext(context);
}



int main() {
    cl_int err;
    const unsigned long IMAGE_WIDTH  = 256*2;
    const unsigned long IMAGE_HEIGHT = 256*2;
    const int SPP = 128*8;//959; //Samples per pixels
    const int ITERATIONS = 8; //Ray recursion count


    const int OBJS_COUNT = 128;
    cl_float16 objs[OBJS_COUNT];
    const int WRAPPED_OBJS_COUNT = 128;
    cl_float16 wrapped_objs[OBJS_COUNT];
    const int LIGHTS_COUNT = 128;
    cl_float16 lights[LIGHTS_COUNT];
    const int MATS_COUNT = 128;
    cl_float16 mats[MATS_COUNT];
    const int TEXS_COUNT = 128;
    cl_float16 texs[TEXS_COUNT];
    int allocated_objects = 0;
    int allocated_lights = 0;


    err = setup_world(&objs[0], &wrapped_objs[0], &lights[0], &mats[0], &texs[0], OBJS_COUNT, WRAPPED_OBJS_COUNT, LIGHTS_COUNT, MATS_COUNT, TEXS_COUNT, &allocated_objects, &allocated_lights);


    scene_settings settings;
    settings.res_width         = 256*2;
    settings.res_height        = settings.res_width;
    settings.samples_per_pixel = 128*8;
    settings.iterations        = 8;


    scene s;
    s.objs_max = 128;
    s.objs = objs;
    s.objs_count = allocated_objects;
    s.wrapped_objs_max = 128;
    s.wrapped_objs = wrapped_objs;
    s.lights_max = 128;
    s.lights = lights;
    s.lights_count = allocated_lights;
    s.mats_max = 128;
    s.mats = mats;
    s.texs_max = 128;
    s.texs = texs;



    render_scene(&s, &settings);


    return 0;
}

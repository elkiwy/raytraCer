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





typedef struct camera {
   float origin[3];
   float lower_left_corner[3];
   float horizontal[3];
   float vertical[3];
} camera;


void init_camera(camera* c) {
   float aspect_ratio = 1.0 / 1.0;
   float viewport_height = 2.0;
   float viewport_width = aspect_ratio * viewport_height;
   float focal_length = 1.0;

   c->origin[0] = 0;
   c->origin[1] = 0;
   c->origin[2] = 0;
   c->horizontal[0] = viewport_width;
   c->horizontal[1] = 0;
   c->horizontal[2] = 0;
   c->vertical[0] = 0.0;
   c->vertical[1] = viewport_height;
   c->vertical[2] = 0.0;

   c->lower_left_corner[0] = c->origin[0] - c->horizontal[0]/2.0 - c->vertical[0]/2.0 - 0;
   c->lower_left_corner[1] = c->origin[1] - c->horizontal[1]/2.0 - c->vertical[1]/2.0 - 0;
   c->lower_left_corner[2] = c->origin[2] - c->horizontal[2]/2.0 - c->vertical[2]/2.0 - focal_length;
}

float random_double_unit() {return rand() / (RAND_MAX + 1.0);}
float random_double(float min, float max) {return min + (max-min)*random_double_unit();}




cl_float3 random_in_unit_sphere(){
    while (1) {
        cl_float3 p = {{random_double(-1,1), random_double(-1,1), random_double(-1,1)}};
        if((p.s[0]*p.s[0] + p.s[1]*p.s[1] + p.s[2]*p.s[2]) < 1){
            return p;
        }
    }
}

#define IMAGE_WIDTH  256*4
#define IMAGE_HEIGHT 256*4

#define CHUNKS_SQRT 2
#define CHUNKS_WIDTH  IMAGE_WIDTH / CHUNKS_SQRT
#define CHUNKS_HEIGHT IMAGE_HEIGHT / CHUNKS_SQRT

#define SPP 128 //Samples per pixels

int main() {
    /**
    ** OpenCL Configuration
     **/
    //Update random seed
    struct timeval t; gettimeofday(&t, NULL); srand(t.tv_usec);

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
    err = clSetKernelArg(kernel, argc, sizeof(cl_uint8), &chunk_data);
    const int ARGC_CHUNK_DATA = argc;
    argc++;

    //Ray Pools = {0,1,2 -> Origin; 3,4,5 -> Direction; 6 -> Used flag; 7 -> ??}
    camera cam; init_camera(&cam);
    const int RAY_POOL_SIZE = (CHUNKS_WIDTH * CHUNKS_HEIGHT) * SPP;
    cl_float8* ray_pool = malloc(sizeof(cl_float8)*RAY_POOL_SIZE);
    cl_mem ray_pool_buffer  = clCreateBuffer(context, F_RW_U, RAY_POOL_SIZE * sizeof(cl_float8), ray_pool, &err);
    if(err < 0) {perror("Couldn't create raypool a buffer"); exit(1);};
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem),   &ray_pool_buffer);
    argc++;

    //Random in unit sphere
    const int RUS_COUNT = 1024;
    cl_float3 rus[RUS_COUNT];
    for (int i=0; i<RUS_COUNT; ++i){rus[i] = random_in_unit_sphere();}
    cl_mem rus_buffer = clCreateBuffer(context, F_R_C, RUS_COUNT * sizeof(cl_float3), rus, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &rus_buffer);
    argc++;

    //Random seeds
    const int RANDOM_SEEDS_COUNT = CHUNKS_HEIGHT*CHUNKS_WIDTH*0.05;
    cl_int random_seeds[RANDOM_SEEDS_COUNT];
    for (int i=0; i<RANDOM_SEEDS_COUNT; ++i){random_seeds[i] = rand();}
    cl_mem random_seeds_buffer = clCreateBuffer(context, F_R_C, RANDOM_SEEDS_COUNT * sizeof(cl_int), random_seeds, &err);
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &random_seeds_buffer);
    argc++;

    //Pass number
    cl_uint pass = 0;
    err = clSetKernelArg(kernel, argc, sizeof(cl_uint),  &pass);
    const int ARGC_PASS = argc;
    argc++;

    //Pixels
    const int PIXEL_COUNT = IMAGE_WIDTH * IMAGE_HEIGHT;
    cl_float4 pixels[PIXEL_COUNT];
    cl_mem pixels_buffer = clCreateBuffer(context, F_RW_C, PIXEL_COUNT * sizeof(cl_float4), pixels, &err);
    if(err < 0) {perror("Couldn't create pixels a buffer"); exit(1);};
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &pixels_buffer);




    /**
    ** Pass execution
     **/
    for (int chunk_y = 0; chunk_y < CHUNKS_SQRT; chunk_y++){
        for (int chunk_x = 0; chunk_x < CHUNKS_SQRT; chunk_x++){
            //Update chunk data
            chunk_data.s[0] = chunk_x;
            chunk_data.s[1] = chunk_y;
            clSetKernelArg(kernel, ARGC_CHUNK_DATA, sizeof(cl_uint8), &chunk_data); // pass

            const float INV_W = 1.0f/(IMAGE_WIDTH-1);
            const float INV_H = 1.0f/(IMAGE_HEIGHT-1);
            for (int py=0; py<CHUNKS_HEIGHT; ++py){
                for (int px=0; px<CHUNKS_WIDTH; ++px){
                    for (int s=0;s<SPP;++s){
                        float u = ((float)(px + CHUNKS_WIDTH*chunk_x) + random_double_unit()) * INV_W;
                        float v = ((float)(py + CHUNKS_HEIGHT*chunk_y) + random_double_unit()) * INV_H;

                        int ray_index = ((py * CHUNKS_WIDTH + px) * SPP) + s;
                        ray_pool[ray_index].s[0] = cam.origin[0];
                        ray_pool[ray_index].s[1] = cam.origin[1];
                        ray_pool[ray_index].s[2] = cam.origin[2];
                        ray_pool[ray_index].s[3] = cam.lower_left_corner[0] + u*cam.horizontal[0] + v*cam.vertical[0] - cam.origin[0];
                        ray_pool[ray_index].s[4] = cam.lower_left_corner[1] + u*cam.horizontal[1] + v*cam.vertical[1] - cam.origin[1];
                        ray_pool[ray_index].s[5] = cam.lower_left_corner[2] + u*cam.horizontal[2] + v*cam.vertical[2] - cam.origin[2];
                        ray_pool[ray_index].s[6] = 0;
                        ray_pool[ray_index].s[7] = 0;
                    }
                }
            }

            for (int i=0; i<16; ++i){
                //Update pass number
                pass = (cl_int)i;
                clSetKernelArg(kernel, ARGC_PASS, sizeof(cl_uint), &pass);

                // Deploy the kernel to a device.
                err = clEnqueueNDRangeKernel(queue, kernel, 2, NULL, global_size, local_size, 0, NULL, NULL);
                //Synchronize en
                err = clEnqueueReadBuffer(queue, pixels_buffer, CL_TRUE, 0, sizeof(pixels), pixels, 0, NULL, NULL);
                printf("Chunk %d %d pass %d done.\n", chunk_x, chunk_y, pass);
            }
            printf("Chunk %d %d done.\n", chunk_x, chunk_y);
        }
    }




















    /**
    ** Kernel output read and processing
     **/
    //Print center pixel
    int c = (IMAGE_WIDTH*(IMAGE_HEIGHT/2)) - IMAGE_WIDTH*0.5;
    printf("Center pixel: %f %f %f\n", pixels[c].s[0], pixels[c].s[1], pixels[c].s[2]);

    //Process the output into an image
    const int CHANNELS = 3;
    unsigned char pixelBytes[PIXEL_COUNT * 3];
    for (int i=0; i<PIXEL_COUNT; i++){
        float samples = pixels[i].s[3];
        pixelBytes[(i*3)+0] = (unsigned char)((pixels[i].s[0]/samples)*255);
        pixelBytes[(i*3)+1] = (unsigned char)((pixels[i].s[1]/samples)*255);
        pixelBytes[(i*3)+2] = (unsigned char)((pixels[i].s[2]/samples)*255);
    }

    //Flip the final image and save it
    stbi_flip_vertically_on_write(1);
    stbi_write_png("test_opencl.png", IMAGE_WIDTH, IMAGE_HEIGHT, CHANNELS, pixelBytes, IMAGE_WIDTH * CHANNELS);


    //Deallocate resources
    free(ray_pool);
    clReleaseKernel(kernel);
    clReleaseMemObject(pixels_buffer);
    clReleaseMemObject(ray_pool_buffer);
    clReleaseMemObject(random_seeds_buffer);
    clReleaseMemObject(rus_buffer);
    clReleaseCommandQueue(queue);
    clReleaseProgram(program);
    clReleaseContext(context);
    return 0;
}

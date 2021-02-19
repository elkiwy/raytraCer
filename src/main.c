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






#define IMAGE_WIDTH  256*4
#define IMAGE_HEIGHT 256*4

#define CHUNKS_SQRT 2
#define CHUNKS_WIDTH  IMAGE_WIDTH / CHUNKS_SQRT
#define CHUNKS_HEIGHT IMAGE_HEIGHT / CHUNKS_SQRT

#define SPP 128 //Samples per pixels
#define ITERATIONS 16 //Ray recursion count

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
    err = clSetKernelArg(kernel, argc, sizeof(cl_uint8), &chunk_data);
    const int ARGC_CHUNK_DATA = argc;
    argc++;

    //Parameters_data = {SPP, RANDOM_COUNT}
    const int RANDOM_SEEDS_COUNT = CHUNKS_HEIGHT*CHUNKS_WIDTH*0.5;
    const int RUS_COUNT = 1024*4;
    cl_int4 parameters_data = {{SPP, RANDOM_SEEDS_COUNT, RUS_COUNT}};
    err = clSetKernelArg(kernel, argc, sizeof(cl_int4), &parameters_data);
    const int ARGC_PARAMETERS_DATA = argc;
    argc++;

    //Ray Pools = {0,1,2 -> Origin; 3,4,5 -> Direction; 6 -> Used flag; 7 -> ??}
    cl_float3 from = {{13, 2, 3}};
    cl_float3   to = {{ 0, 0, 0}};
    cl_float3  vup = {{ 0, 1, 0}};
    float dist_to_focus = 10.0;
    float aperture = 0.1;
    camera cam; init_camera(&cam, from, to, vup, 20.0, 1.0, aperture, dist_to_focus);
    const int RAY_POOL_SIZE = (CHUNKS_WIDTH * CHUNKS_HEIGHT) * SPP;
    cl_float16* ray_pool = malloc(sizeof(cl_float16)*RAY_POOL_SIZE);
    cl_mem ray_pool_buffer  = clCreateBuffer(context, F_RW_U, RAY_POOL_SIZE * sizeof(cl_float16), ray_pool, &err);
    if(err < 0) {perror("Couldn't create raypool a buffer"); exit(1);};
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem),   &ray_pool_buffer);
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

    //Iterations data
    cl_uint2 iteration_data = (cl_uint2){{0, ITERATIONS}};
    err = clSetKernelArg(kernel, argc, sizeof(cl_uint2),  &iteration_data);
    const int ARGC_ITERATIONS_DATA = argc;
    argc++;

    //Pixels
    const int PIXEL_COUNT = IMAGE_WIDTH * IMAGE_HEIGHT;
    const int OUTPUT_SIZE = PIXEL_COUNT;
    cl_float4 output[PIXEL_COUNT];
    cl_mem output_buffer = clCreateBuffer(context, F_RW_C, OUTPUT_SIZE * sizeof(cl_float4), output, &err);
    if(err < 0) {perror("Couldn't create pixels a buffer"); exit(1);};
    err = clSetKernelArg(kernel, argc, sizeof(cl_mem), &output_buffer);


    printf("Kernel ready for execution:\n");fflush(stdout);


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
                        ray_pool[ray_index] = get_ray(&cam, u, v);
                    }
                }
            }

            //Ensure the buffer gets updated with the new set of ray each iteration (not necessary on certain GPU, but needed for others)
            err = clEnqueueWriteBuffer(queue, ray_pool_buffer, CL_FALSE, 0, RAY_POOL_SIZE * sizeof(cl_float16), ray_pool, 0, NULL, NULL);

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
            printf("Chunk %d %d done.\n", chunk_x, chunk_y);
        }
    }




















    /**
    ** Kernel output read and processing
     **/
    //Process the output into an image
    const int CHANNELS = 3;
    unsigned char pixelBytes[PIXEL_COUNT * 3];
    for (int i=0; i<PIXEL_COUNT; i++){


        float samples_for_iteration = output[i].s[3];
        float r = output[i].s[0]/samples_for_iteration;
        float g = output[i].s[1]/samples_for_iteration;
        float b = output[i].s[2]/samples_for_iteration;

        if (i==0  || i == 201472 || i == 54010 || i==PIXEL_COUNT-1){printf("\n\n%d pixel iteration %d: %f %f %f (%f samples)\n", i, 0, output[i].s[0], output[i].s[1], output[i].s[2], output[i].s[3]);}
        if (i==0  || i == 201472 || i == 54010 || i==PIXEL_COUNT-1){printf("%d final values : %f %f %f \n", i, r, g, b);}

        r = sqrt(r);
        g = sqrt(g);
        b = sqrt(b);

        if (i==0 || i == 201472 || i == 54010 || i==PIXEL_COUNT-1){printf("%d final values after gamma: %f %f %f \n", i, r, g, b);}


        //float samples = output[i].s[3];
        //float scale = 1.0f/samples;
        pixelBytes[(i*3)+0] = (unsigned char)(r*255);
        pixelBytes[(i*3)+1] = (unsigned char)(g*255);
        pixelBytes[(i*3)+2] = (unsigned char)(b*255);
    }

    //Flip the final image and save it
    stbi_flip_vertically_on_write(1);
    stbi_write_png("test_opencl.png", IMAGE_WIDTH, IMAGE_HEIGHT, CHANNELS, pixelBytes, IMAGE_WIDTH * CHANNELS);


    //Deallocate resources
    free(ray_pool);
    clReleaseKernel(kernel);
    clReleaseMemObject(output_buffer);
    clReleaseMemObject(ray_pool_buffer);
    clReleaseMemObject(random_seeds_buffer);
    clReleaseMemObject(rus_buffer);
    clReleaseCommandQueue(queue);
    clReleaseProgram(program);
    clReleaseContext(context);
    return 0;
}

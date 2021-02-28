#include "utils.h"

#include "kernelSource.h"


/**
 * Randoms
 * */

float random_double_unit() {
    return rand() / (RAND_MAX + 1.0);
}



float random_double(float min, float max) {
    return min + (max-min)*random_double_unit();
}



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





/**
 * Vector Math
 * */

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








/**
 * Perlin
 * */

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










/**
 * OpenCL
 * */

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
cl_program build_program(cl_context ctx, cl_device_id dev) {
    //Create program from file
    int err;
    char* program_buffer = src_program_cl;
    size_t program_size = src_program_cl_len;
    cl_program program = clCreateProgramWithSource(ctx, 1, (const char**)&program_buffer, &program_size, &err);
    if(err < 0) {perror("Couldn't create the program"); exit(1);}

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









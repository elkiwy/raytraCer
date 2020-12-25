#include <stdio.h>
#include <math.h>



/**
 *
 * Vec3
 *
 * */

//Struct
typedef struct{double x; double y; double z;}vec3;
inline static void vec3_print(vec3* v){printf("Vec3: %f %f %f\n", v->x, v->y, v->z);}

//Vector operations without copies
inline static vec3 vec3_sum(vec3* v1, vec3* v2){vec3 result = {v1->x + v2->x, v1->y + v2->y, v1->z + v2->z}; return result;}
inline static vec3 vec3_sub(vec3* v1, vec3* v2){vec3 result = {v1->x - v2->x, v1->y - v2->y, v1->z - v2->z}; return result;}
inline static vec3 vec3_mul(vec3* v1, vec3* v2){vec3 result = {v1->x * v2->x, v1->y * v2->y, v1->z * v2->z}; return result;}
inline static vec3 vec3_mul_k(vec3* v, double k){vec3 result = {v->x * k, v->y * k, v->z * k}; return result;}
inline static vec3 vec3_div_k(vec3* v, double k){vec3 result = {v->x * 1/k, v->y * 1/k, v->z * 1/k}; return result;}
inline static double vec3_length(vec3* v){return sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));}
inline static double vec3_length_squared(vec3* v){return (v->x * v->x) + (v->y * v->y) + (v->z * v->z);}
inline static double vec3_dot(vec3* v1, vec3* v2){return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;}
inline static vec3 vec3_cross(vec3* v1, vec3* v2){vec3 result = {v1->y * v2->z - v1->z * v2->y, v1->z * v2->x - v1->x * v2->z, v1->x * v2->y - v1->y * v2->x}; return result;}
inline static vec3 vec3_unit(vec3* v1){vec3 res = vec3_div_k(v1, vec3_length(v1));return res;}

//Vector operation with copies
inline static vec3 vec3c_sum(vec3 v1, vec3 v2){vec3 result = {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z}; return result;}
inline static vec3 vec3c_sub(vec3 v1, vec3 v2){vec3 result = {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z}; return result;}
inline static vec3 vec3c_mul(vec3 v1, vec3 v2){vec3 result = {v1.x * v2.x, v1.y * v2.y, v1.z * v2.z}; return result;}
inline static vec3 vec3c_mul_k(vec3 v, double k){vec3 result = {v.x * k, v.y * k, v.z * k}; return result;}
inline static vec3 vec3c_div_k(vec3 v, double k){vec3 result = {v.x * 1/k, v.y * 1/k, v.z * 1/k}; return result;}
inline static double vec3c_length(vec3 v){return sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));}
inline static double vec3c_length_squared(vec3 v){return (v.x * v.x) + (v.y * v.y) + (v.z * v.z);}
inline static double vec3c_dot(vec3 v1, vec3 v2){return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;}
inline static vec3 vec3c_cross(vec3 v1, vec3 v2){vec3 result = {v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x}; return result;}
inline static vec3 vec3c_unit(vec3 v1){vec3 res = vec3c_div_k(v1, vec3c_length(v1));return res;}

//Names
typedef vec3 point3;
typedef vec3 color;




/**
 *
 * Ray
 *
 * */

typedef struct{point3 orig; vec3 dir;} ray;

inline static point3 ray_at(ray* r, double t){
    vec3 a = vec3_mul_k(&r->dir, t);
    return vec3_sum(&r->orig, &a);
}












/**
 *
 * Color
 *
 * */


void write_color(FILE* f, color* pixel_color){
    fprintf(f, "%d %d %d\n", (int)(255.999*pixel_color->x), (int)(255.999*pixel_color->y), (int)(255.999*pixel_color->z));
}










/**
 *
 * Main
 *
 * */

double hit_sphere(point3* center, double radius, ray* r){
    vec3 oc = vec3_sub(&(r->orig), center);
    double a = vec3_dot(&(r->dir), &(r->dir));
    double b = 2.0 * vec3_dot(&oc, &(r->dir));
    double c = vec3_dot(&oc, &oc) - radius * radius;
    double discriminant = b*b - 4*a*c;

    if (discriminant < 0){
        return -1.0;
    }else{
        return (-b - sqrt(discriminant)) / (2.0*a);
    }
}

color ray_color(ray* r){
    point3 sphere_center = {0,0,-1};
    double t = hit_sphere(&sphere_center, 0.5, r);
    if (t>0){
        vec3 N = vec3c_unit(vec3c_sub(ray_at(r, t), sphere_center));
        color c = {N.x + 1, N.y + 1, N.z + 1};
        return vec3_mul_k(&c, 0.5);
    }

    vec3 unit_dir = vec3_unit(&r->dir);
    t = 0.5 * (unit_dir.y + 1.0);
    color white = {1.0, 1.0, 1.0};
    color cyan = {0.5, 0.7, 1.0};
    return vec3c_sum(vec3_mul_k(&white, (1.0 - t)), vec3_mul_k(&cyan, t));
}











int main(int argc, char** argv) {
    //Output file
    char* outputFilePath = argv[1];
    char buffer[1024*1024]; // 1 MB buffer
    FILE* f = fopen(outputFilePath, "w+");

    //Define camera
    const double ASPECT_RATIO = 16.0/9.0;
    const int IMAGE_WIDTH = 400;
    const int IMAGE_HEIGHT = (int)(IMAGE_WIDTH / ASPECT_RATIO);

    double viewport_height = 2.0;
    double viewport_width = ASPECT_RATIO * viewport_height;
    double focal_length = 1.0;

    point3 origin = {0, 0, 0};
    vec3 horizontal = {viewport_width, 0, 0};
    vec3 vertical = {0, viewport_height, 0};
    vec3 depth = {0, 0, focal_length};
    vec3 lower_left_corner = vec3c_sub(vec3c_sub(vec3c_sub(origin, vec3_div_k(&horizontal, 2.0)), vec3_div_k(&vertical, 2.0)), depth);



    //Processing
    int progress = 0;
    fprintf(f, "P3\n%d %d\n255\n", IMAGE_WIDTH, IMAGE_HEIGHT);
    for(int j=IMAGE_HEIGHT-1; j>=0; --j){
        //Print progress
        if (progress < (1-((double)j/(double)(IMAGE_HEIGHT-1)))*100){printf("Progress: %d%%\n",progress); progress += 10;}

        for(int i=0; i<IMAGE_WIDTH; ++i){

            double u = (double)i / (double)(IMAGE_WIDTH - 1);
            double v = (double)j / (double)(IMAGE_HEIGHT - 1);

            vec3 dir = vec3c_sub(vec3c_sum(vec3c_sum(lower_left_corner, vec3_mul_k(&horizontal, u)), vec3_mul_k(&vertical, v)), origin);
            ray r = {origin, dir};
            color pixel_color = ray_color(&r);
            write_color(f, &pixel_color);
        }
    }

    printf("Progress: 100%%\n");
    fclose(f);


    return 0;
}

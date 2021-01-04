#include "aabb.h"



aabb* aabb_init(point3 a, point3 b){
    aabb* bb = malloc(sizeof(aabb));
    bb->minimum = a;
    bb->maximum = b;
    return bb;
}

void aabb_free(aabb* b){
    free(b);
}

void aabb_print(aabb* b){
    printf("Bounding box:\n");
    printf("  Min:");
    vec3_print(&b->minimum);
    printf("  Max:");
    vec3_print(&b->maximum);
    fflush(stdout);
}


void swapDoubles(double* a, double*b){double t = *a; *a = *b; *b = t;}

int aabb_hit(aabb* box, ray* r, double tmin, double tmax){
    //Test for X
    float invD = 1.0f/(float)r->dir.x;
    double t0x = (box->minimum.x - r->orig.x)*invD;
    double t1x = (box->maximum.x - r->orig.x)*invD;
    if (invD < 0.0f){swapDoubles(&t0x, &t1x);}
    tmin = t0x > tmin ? t0x : tmin;
    tmax = t1x < tmax ? t1x : tmax;
    if (tmax <= tmin) {return 0;}

    //Test for Y
    invD = 1.0f/(float)r->dir.y;
    double t0y = (box->minimum.y - r->orig.y)*invD;
    double t1y = (box->maximum.y - r->orig.y)*invD;
    if (invD < 0.0f){swapDoubles(&t0y, &t1y);}
    tmin = t0y > tmin ? t0y : tmin;
    tmax = t1y < tmax ? t1y : tmax;
    if (tmax <= tmin) {return 0;}

    //Test for Z
    invD = 1.0f/(float)r->dir.z;
    double t0z = (box->minimum.z - r->orig.z)*invD;
    double t1z = (box->maximum.z - r->orig.z)*invD;
    if (invD < 0.0f){swapDoubles(&t0z, &t1z);}
    tmin = t0z > tmin ? t0z : tmin;
    tmax = t1z < tmax ? t1z : tmax;
    if (tmax <= tmin) {return 0;}

    //Hit success
    return 1;
}


aabb surrounding_box(aabb box0, aabb box1){
    point3 small = {
        fmin(box0.minimum.x, box1.minimum.x),
        fmin(box0.minimum.y, box1.minimum.y),
        fmin(box0.minimum.z, box1.minimum.z)
    };

    point3 big = {
        fmax(box0.maximum.x, box1.maximum.x),
        fmax(box0.maximum.y, box1.maximum.y),
        fmax(box0.maximum.z, box1.maximum.z)
    };

    return (aabb){small, big};
}

#include "ray.h"


///Move the origin point through the ray for T distance
inline point3 ray_at(ray* r, double t){
    vec3 a = vec3_mul_k(&r->dir, t);
    return vec3_sum(&r->orig, &a);
}

///Print a ray object
void ray_print(char* s, ray* r){
    printf("%s: orig: [%.2f %.2f %.2f] dir: [%.2f %.2f %.2f]\n", s, r->orig.x, r->orig.y, r->orig.z, r->dir.x, r->dir.y, r->dir.z);
}

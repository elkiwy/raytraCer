#include "onb.h"




///Initialize a new orthogonal base from a vector w
void onb_init_from_w(onb* o, vec3* n){
    o->axis[2] = vec3_unit(n);
    vec3 a = (fabs(o->axis[2].x) > 0.9) ? (vec3){0,1,0} : (vec3){1,0,0};
    o->axis[1] = vec3c_unit(vec3c_cross(o->axis[2], a));
    o->axis[0] = vec3c_cross(o->axis[2], o->axis[1]);
}


///Get the relative coordinates of a vector in this ONB space
vec3 onb_local(onb* o, double a, double b, double c){
    vec3 va = vec3_mul_k(&o->axis[0], a);
    vec3 vb = vec3_mul_k(&o->axis[1], b);
    vec3 vc = vec3_mul_k(&o->axis[2], c);
    return vec3c_sum(va, vec3c_sum(vb, vc));
}


///Get the relative coordinates of a vector in this ONB space
vec3 onb_localv(onb* o, vec3 a){
    vec3 va = vec3_mul_k(&o->axis[0], a.x);
    vec3 vb = vec3_mul_k(&o->axis[1], a.y);
    vec3 vc = vec3_mul_k(&o->axis[2], a.z);
    return vec3c_sum(va, vec3c_sum(vb, vc));
}

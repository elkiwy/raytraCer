#ifndef __ONB_H_
#define __ONB_H_

#include "utils.h"
#include "vec3.h"

//Struct
typedef struct onb {vec3 axis[3];}onb;

//Init
void onb_init_from_w(onb* o, vec3* n);

//Features
vec3 onb_local(onb* o, double a, double b, double c);
vec3 onb_localv(onb* o, vec3 a);


#endif // __ONB_H_

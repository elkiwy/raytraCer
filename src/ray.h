#ifndef __RAY_H_
#define __RAY_H_


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "vec3.h"

struct hittable_list;


/**
 *
 * Ray
 *
 * */

///Ray structure
typedef struct{point3 orig; vec3 dir; double time;} ray;

///Move the origin point through the ray for T distance
point3 ray_at(ray* r, double t);

///Print a ray object
void ray_print(char* s, ray* r);




///Prints a color in hex RGB
void color_printHex(char* s, color* c);

///Write a pixel color to the output file, the color is averaged between the samples_per_pixel
void write_color(unsigned char* pixels, color* pixel_color, int samples_per_pixel, int index);

///Cast a ray into the world and retrieve the color of that ray
color ray_color(ray* r, color* background, struct hittable_list* world, struct hittable_list* lights, int recursion_depth);

#endif // __RAY_H_

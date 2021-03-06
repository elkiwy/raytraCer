#ifndef __UTILS_H_
#define __UTILS_H_


#include <stdlib.h>


/**
 *
 * Utilities and constants
 *
 * */


// Constants
const double PI = 3.1415926535897932385;

// Utility Functions
inline static double deg2rad(double deg) {return deg * PI / 180.0;}

inline static double random_double(){return rand() / (RAND_MAX + 1.0);}
inline static double random_double_scaled(double min, double max){return min + (max-min)*random_double();}

inline static double clamp(double x, double min, double max){
    return (x < min) ? min : ((x > max) ? max : x);
}



#endif // __UTILS_H_

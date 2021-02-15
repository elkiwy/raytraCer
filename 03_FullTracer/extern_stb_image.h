#ifndef __EXTERN_STB_IMAGE_H_
#define __EXTERN_STB_IMAGE_H_


// Disable pedantic warnings for this external library.
#ifdef _MSC_VER
    // Microsoft Visual C++ Compiler
    #pragma warning (push, 0)
#endif



#define STB_IMAGE_IMPLEMENTATION
#include "external/stb_image.h"




// Restore warning levels.
#ifdef _MSC_VER
    // Microsoft Visual C++ Compiler
    #pragma warning (pop)
#endif



#endif // __EXTERN_STB_IMAGE_H_
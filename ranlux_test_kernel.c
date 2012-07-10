
%(defines)s

#define RANLUXCL_LUX %(luxuaryFactor)d


#ifdef ENABLE_DOUBLE
    #define FLOAT_TYPE double
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #define RANLUXCL_SUPPORT_DOUBLE
#else
    #define FLOAT_TYPE float
#endif

#include "ranlux.cl"


__kernel void ranlux_test_kernel(__global FLOAT_TYPE *randomNumbers)
{
    
}

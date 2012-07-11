
%(defines)s

#define RANLUXCL_LUX %(luxuaryFactor)d


#ifdef ENABLE_DOUBLE
    #define FLOAT_TYPE double
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #define RANLUXCL_SUPPORT_DOUBLE
#else
    #define FLOAT_TYPE float
#endif

#include "pyopencl-ranluxcl.cl"
__kernel void Kernel_Ranluxcl_Init(__private uint ins,
__global ranluxcl_state_t *ranluxcltab)
{
    ranluxcl_initialization(ins, ranluxcltab);
}

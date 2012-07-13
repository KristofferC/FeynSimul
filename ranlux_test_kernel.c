
%(defines)s

#define RANLUXCL_LUX %(luxuaryFactor)s


#ifdef ENABLE_DOUBLE
    #define FLOAT_TYPE double
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #define RANLUXCL_SUPPORT_DOUBLE
#else
    #define FLOAT_TYPE float
#endif

#include "pyopencl-ranluxcl.cl"

__kernel void ranlux_init_kernel(__private uint ins,
                                 __global ranluxcl_state_t *ranluxcltab)
{
    ranluxcl_initialization(ins, ranluxcltab);
}

__kernel void ranlux_test_kernel(__global FLOAT_TYPE *randomsOut,
                                 __global ranluxcl_state_t *ranluxcltab)
{
    uint threadId = get_global_id(0) + get_global_id(1) * get_global_size(0);
    
    //ranluxclstate stores the state of the generator.
    ranluxcl_state_t ranluxclstate;

    //Download state into ranluxclstate struct.
    ranluxcl_download_seed(&ranluxclstate, ranluxcltab);

    //Generate a float4 with each component on (0,1),
    //end points not included. We can call ranluxcl as many
    //times as we like until we upload the state again.
    float4 randomnr;
    uint randOffset;
    uint randLocalOffset = threadId * %(randsPerThread)s;
    
    for (int i = 0; i <= %(randsPerThread)s; i++)
    {
        randOffset = randLocalOffset + 4 * i;
        randomnr = ranluxcl32(&ranluxclstate);
        randomsOut[randOffset + 0] = randomnr.x;
        randomsOut[randOffset + 1] = randomnr.y;
        randomsOut[randOffset + 2] = randomnr.z;
        randomsOut[randOffset + 3] = randomnr.w;
    }
    
    //Upload state again so that we don't get the same
    //numbers over again the next time we use ranluxcl.
    ranluxcl_upload_seed(&ranluxclstate, ranluxcltab);
}

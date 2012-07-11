
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

__kernel void ranlux_test_kernel(__global ranluxcl_state_t *ranluxcltab)
{
    //ranluxclstate stores the state of the generator.
    ranluxcl_state_t ranluxclstate;

    //Download state into ranluxclstate struct.
    ranluxcl_download_seed(&ranluxclstate, ranluxcltab);

    //Generate a float4 with each component on (0,1),
    //end points not included. We can call ranluxcl as many
    //times as we like until we upload the state again.
    float4 randomnr = ranluxcl32(&ranluxclstate);

    //Upload state again so that we don't get the same
    //numbers over again the next time we use ranluxcl.
    ranluxcl_upload_seed(&ranluxclstate, ranluxcltab);
}

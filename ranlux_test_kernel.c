%(defines)s

#ifdef ENABLE_DOUBLE
    #define DATA_TYPE double
    #define DATA_TYPE_V double4
    #define DATA_TYPE_F ranluxcl64
    #define FLOAT_TYPE double
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    
    #ifdef ENABLE_RANLUX
        #define RANLUXCL_SUPPORT_DOUBLE
    #endif
#else
    #ifdef USE_RANLUXCL_INT
        #define DATA_TYPE uint
        #define DATA_TYPE_V uint4
        #define DATA_TYPE_F ranluxclint
    #else
        #define DATA_TYPE float
        #define DATA_TYPE_V float4
        #define DATA_TYPE_F ranluxcl32
    #endif
    
    #define FLOAT_TYPE float
#endif


//__kernel void ranlux_init_kernel(__global uint *ins,
//                                 __global ranluxcl_state_t *ranluxcltab)
//{
//    ranluxcl_initialization(ins, ranluxcltab);
//}


#ifdef ENABLE_RANLUX
uint4 ranluxclint(ranluxcl_state_t *ranluxclstate)
{
    return convert_uint4(ranluxcl32(ranluxclstate) * (float4)1.6777216E7f);
}

#define RANLUXCL_LUX %(luxuaryFactor)s
#include "pyopencl-ranluxcl.cl"
__kernel void ranlux_test_kernel(__global uint *ins,
#ifdef RETURN_RANDOMS
                                 __global DATA_TYPE *randomsOut,
#endif
                                 __global ranluxcl_state_t *ranluxcltab)
{
    ranluxcl_initialization(ins, ranluxcltab);
    
    uint threadId = get_global_id(0) + get_global_id(1) * get_global_size(0);
    
    //ranluxclstate stores the state of the generator.
    ranluxcl_state_t ranluxclstate;

    //Download state into ranluxclstate struct.
    ranluxcl_download_seed(&ranluxclstate, ranluxcltab);

    //Generate a float4 with each component on (0,1),
    //end points not included. We can call ranluxcl as many
    //times as we like until we upload the state again.'
    DATA_TYPE_V randomnr;
    
    uint randOffset;
    uint randLocalOffset = threadId * %(randsPerThread)s;
    
    for (int i = 0; i <= %(randsPerThread)s / 4; i++)
    {
        randOffset = randLocalOffset + 4 * i;
        randomnr = DATA_TYPE_F (&ranluxclstate);

#ifdef RETURN_RANDOMS
        randomsOut[randOffset + 0] = randomnr.x;
        randomsOut[randOffset + 1] = randomnr.y;
        randomsOut[randOffset + 2] = randomnr.z;
        randomsOut[randOffset + 3] = randomnr.w;
#endif
    }
    
    //Upload state again so that we don't get the same
    //numbers over again the next time we use ranluxcl.
    ranluxcl_upload_seed(&ranluxclstate, ranluxcltab);
}
#endif


#ifndef ENABLE_RANLUX
inline void xorshift (uint4 *seedPtr)
{
    uint t;
    t = (*seedPtr).x ^ ((*seedPtr).x << 11);
    (*seedPtr).xyz=(*seedPtr).yzw;
    (*seedPtr).w = ((*seedPtr).w ^ ((*seedPtr).w >> 19) ^ (t ^ (t >> 8)));
}

inline FLOAT_TYPE randFloat(uint4 *seedPtr)
{
    xorshift(seedPtr);
    return (*seedPtr).w * 2.328306437080797e-10;
    //                    5.421010862427522e-10
    //return ((double) (*seedPtr).w) * 2.32830643653869628906e-10;
}

__kernel void xorshift_test_kernel(__global uint *ins
#ifdef RETURN_RANDOMS
                                 , __global FLOAT_TYPE *randomsOut
#endif
                                  )
{
    uint threadId = get_global_id(0) + get_global_id(1) * get_global_size(0);
    
    uint4 seed,seedG;
    seed.x = ins[threadId * 4 + 0];
    seed.y = ins[threadId * 4 + 1];
    seed.z = ins[threadId * 4 + 2];
    seed.w = ins[threadId * 4 + 3];

    uint randOffset = threadId * %(randsPerThread)s;
    
    FLOAT_TYPE randomnr;
    
    for (int i = 0; i <= %(randsPerThread)s; i++)
    {
        randomnr = randFloat(&seed);
#ifdef RETURN_RANDOMS
        randomsOut[randOffset + i] = randomnr;
#endif
    }
    
    ins[threadId * 4] = seed.x;
    ins[threadId * 4 + 1] = seed.y;
    ins[threadId * 4 + 2] = seed.z;
    ins[threadId * 4 + 3] = seed.w;
}
#endif

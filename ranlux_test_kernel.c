%(defines)s

#ifdef ENABLE_DOUBLE
    #define FLOAT_TYPE double
    
    #ifdef ENABLE_RANLUX
        #pragma OPENCL EXTENSION cl_khr_fp64 : enable
        #define RANLUXCL_SUPPORT_DOUBLE
    #endif
#else
    #define FLOAT_TYPE float
#endif

#ifdef ENABLE_RANLUX
    #define RANLUXCL_LUX %(luxuaryFactor)s
    #include "pyopencl-ranluxcl.cl"
#endif

//__kernel void ranlux_init_kernel(__global uint *ins,
//                                 __global ranluxcl_state_t *ranluxcltab)
//{
//    ranluxcl_initialization(ins, ranluxcltab);
//}

inline void xorshift (uint4 *seedPtr)
{
    uint t;
    t = (*seedPtr).x ^ ((*seedPtr).x << 11);
    (*seedPtr).xyz=(*seedPtr).yzw;
    (*seedPtr).w = ((*seedPtr).w ^ ((*seedPtr).w >> 19) ^ (t ^ (t >> 8)));
}

inline FLOAT_TYPE
randFloat(uint4 *seedPtr)
{
    xorshift(seedPtr);
    return (*seedPtr).w * 2.328306437080797e-10;
}

__kernel void ranlux_test_kernel(__global uint *ins,
#ifdef RETURN_RANDOMS
                                 __global FLOAT_TYPE *randomsOut,
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
#ifdef ENABLE_DOUBLE
    double4 randomnr;
#else
    float4 randomnr;
#endif
    uint randOffset;
    uint randLocalOffset = threadId * %(randsPerThread)s;
    
    for (int i = 0; i <= %(randsPerThread)s / 4; i++)
    {
        randOffset = randLocalOffset + 4 * i;
#ifdef ENABLE_DOUBLE
        randomnr = ranluxcl64(&ranluxclstate);
#else
        randomnr = ranluxcl32(&ranluxclstate);
#endif

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

__kernel void xorshift_test_kernel(__global uint *ins,
#ifdef RETURN_RANDOMS
                                 __global FLOAT_TYPE *randomsOut,
#endif
                                  )
{
    uint threadId = get_global_id(0) + get_global_id(1) * get_global_size(0);
    
    uint4 seed,seedG;
    seed.x = ins[threadId * 4 + 0];
    seed.y = ins[threadId * 4 + 1];
    seed.z = ins[threadId * 4 + 2];
    seed.w = ins[threadId * 4 + 3];

    uint randLocalOffset = threadId * %(randsPerThread)s;
    
    FLOAT_TYPE randomnr;
    
    for (int i = 0; i <= %(randsPerThread)s; i++)
    {
        randomnr = randFloat (&seed);
#ifdef RETURN_RANDOMS
        randomsOut[randOffset + i] = randomnr;
#endif
    }
    
    seeds[threadId * 4] = seed.x;
    seeds[threadId * 4 + 1] = seed.y;
    seeds[threadId * 4 + 2] = seed.z;
    seeds[threadId * 4 + 3] = seed.w;
}

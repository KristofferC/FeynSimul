%(defines)s

#ifdef ENABLE_DOUBLE
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    
    #define DATA_TYPE double
    #define DATA_TYPE_V double4
    #define DATA_TYPE_F ranluxcl64(ranluxclstate)
    #define FLOAT_TYPE double
    #define FLOAT_TYPE_V double4
    #define RANLUX_WRAPPER ranluxWrapper(&ranluxclstate, &randCountL, &random_temp, this_run_i)
    
    #ifdef ENABLE_RANLUX
        #define RANLUXCL_SUPPORT_DOUBLE
    #endif
#else
    #ifdef USE_RANLUXCL_INT
        #define DATA_TYPE uint
        #define DATA_TYPE_V uint4
        #define RANLUX_WRAPPER ranluxIntWrapper(&ranluxclstate, &randCountL, &random_temp, this_run_i)
    #else
        #define DATA_TYPE float
        #define DATA_TYPE_V float4
        #define RANLUX_WRAPPER ranluxWrapper(&ranluxclstate, &randCountL, &random_temp, this_run_i)
    #endif
    
    #define DATA_TYPE_F ranluxcl32(ranluxclstate)
    #define FLOAT_TYPE float
    #define FLOAT_TYPE_V float4
#endif




#ifdef ENABLE_RANLUX
#define RANLUXCL_LUX %(luxuaryFactor)s
#include "pyopencl-ranluxcl.cl"

inline uint4 ranluxclint(ranluxcl_state_t *ranluxclstate)
{
    return convert_uint4(ranluxcl32(ranluxclstate) * (float4) %(ranluxIntMax)s);
}

//__kernel void ranlux_init_kernel(__global uint *ins,
//                                 __global ranluxcl_state_t *ranluxcltab)
//{
//    ranluxcl_initialization(ins, ranluxcltab);
//}
    
__kernel void ranlux_test_kernel_init(__global uint *ins,
                                      __global ranluxcl_state_t *ranluxcltab)
{
    // global group id: get_group_id(0) + get_global_id(1) * get_global_size(0) / get_local_size(0)
    ranluxcl_init((ulong) *(ins + get_group_id(0) + (get_global_id(1) * get_global_size(0)) / get_local_size(0)), ranluxcltab + (get_global_id(0) + get_global_id(1) * get_global_size(0) + get_global_id(2) * get_global_size(0) * get_global_size(1)));
    //ranluxcl_initialization(*ins, ranluxcltab);
}

union ranlux_vector_union
{
    FLOAT_TYPE_V ranlux_vector;
    FLOAT_TYPE ranlux_array[4];
};

inline FLOAT_TYPE ranluxWrapper(ranluxcl_state_t *ranluxclstate,
                                uint *randCount,
                                union ranlux_vector_union *random_temp,
                                uint this_run_i)
{
    if(((*randCount) & 3) == 0 || this_run_i == 0)
    {
        (*random_temp).ranlux_vector = DATA_TYPE_F;// (ranluxclstate);
    }
    
    return (*random_temp).ranlux_array[((*randCount)++) & 3];
}

inline uint ranluxIntWrapper(ranluxcl_state_t *ranluxclstate,
                            uint *randCount,
                            union ranlux_vector_union *random_temp,
                            uint this_run_i)
{
    if(((*randCount) & 3) == 0 || this_run_i == 0)
    {
        (*random_temp).ranlux_vector = DATA_TYPE_F;// (ranluxclstate);
    }
    
    return (uint) ((*random_temp).ranlux_array[((*randCount)++) & 3] * ((FLOAT_TYPE) %(ranluxIntMax)s));
}

__kernel void ranlux_test_kernel(__global uint *ins,
#ifdef RETURN_RANDOMS
                                 __global DATA_TYPE *randomsOut,
#endif
                                 __global uint *testOut,
                                 __global ranluxcl_state_t *ranluxcltab,
                                 __global uint *randCountG)
{
    uint threadId = get_global_id(0) + get_global_id(1) * get_global_size(0);
    
    //ranluxclstate stores the state of the generator.
    ranluxcl_state_t ranluxclstate;

    //ranluxcl_initialization(*ins, ranluxcltab);

    //Download state into ranluxclstate struct.
    ranluxcl_download_seed(&ranluxclstate, ranluxcltab);
    __private uint randCountL = randCountG[threadId];
    //*(&ranluxclstate) = ranluxcltab[get_global_id(0) + get_global_id(1) * get_global_size(0) + get_global_id(2) * get_global_size(0) * get_global_size(1)];

    //Generate a float4 with each component on (0,1),
    //end points not included. We can call ranluxcl as many
    //times as we like until we upload the state again.'

    DATA_TYPE randomnr;

    union ranlux_vector_union random_temp;

    uint randOffset;
    uint randLocalOffset = threadId * %(randsPerThread)s;
    uint this_run_i;
    
    for (int i = 0; i < %(randsPerThread)s; i++)
    {
        this_run_i = i;
        randOffset = randLocalOffset + 4 * i;
        randomnr = RANLUX_WRAPPER; // (&ranluxclstate, &randCount, &random_temp);
        //randomnr = DATA_TYPE_F (&ranluxclstate);
        testOut[randLocalOffset + i] = get_group_id(0) + get_global_id(1) * get_global_size(0) / get_local_size(0);
        
#ifdef RETURN_RANDOMS
        //randomsOut[randOffset + 0] = randomnr.x;
        //randomsOut[randOffset + 1] = randomnr.y;
        //randomsOut[randOffset + 2] = randomnr.z;
        //randomsOut[randOffset + 3] = randomnr.w;
        randomsOut[randLocalOffset + i] = randomnr;
        //randomsOut[randLocalOffset + i] = (get_global_id(0) + get_global_id(1) * get_global_size(0) + get_global_id(2) * get_global_size(0) * get_global_size(1));
        //randomsOut[randLocalOffset + i] = *(ins + get_group_id(0) + get_global_id(1) * get_global_size(0) / get_local_size(0));
        #endif
    /*
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
#endif*/
    }
    
    //Upload state again so that we don't get the same
    //numbers over again the next time we use ranluxcl.
    ranluxcl_upload_seed(&ranluxclstate, ranluxcltab);
    //ranluxcltab[get_global_id(0) + get_global_id(1) * get_global_size(0) + get_global_id(2) * get_global_size(0) * get_global_size(1)] = *(&ranluxclstate);
    randCountG[threadId] = randCountL;
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

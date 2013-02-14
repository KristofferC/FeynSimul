// This file is part of FeynSimul.
//
// FeynSimul is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// FeynSimul is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FeynSimul.  If not, see <http://www.gnu.org/licenses/>.

//##############################################################################
//kernel.c
//______________________________________________________________________________
//Imported by: ../host.py
//______________________________________________________________________________
//Description: This is the OpenCL kernel and is the code that runs on the GPU.
//             This contains the Monte Carlo calculations that are to be per-
//             performed and some helper methods for generating random numbers
//             numbers with the Xorshift PRNG.
//______________________________________________________________________________
//Authors:     Kristoffer Carlsson, Patric Holmvall, Petter Saterskog
//Date:        2011-12-25
//Licensing:   GPL
//##############################################################################

#define eps %(epsilon)s

//##############################################################################
//#                                 Defines                                    #
//##############################################################################
//Description: This is a placeholder for defines used in the program.
#ifdef ENABLE_DOUBLE
    // Set float precision to double
    #define FLOAT_TYPE double
    #define FLOAT_TYPE_VECTOR double4
    
    // Check that pragmas for 64bit actually exists
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    //#ifdef cl_khr_fp64
    //    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    //#elif defined(cl_amd_fp64)
    //    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    //#else
    //    #error "Double precision floating point not supported by OpenCL implementation."
    //#endif
#else
    // Set float precision to single
    #define FLOAT_TYPE float
    #define FLOAT_TYPE_VECTOR float4
#endif

%(defines)s

#ifdef ENABLE_GLOBAL_PATH
	#define PATH_TYPE_KEYWORD __global
#else
	#ifdef ENABLE_PARALELLIZE_PATH
		#define PATH_TYPE_KEYWORD __local
	#else
		#define PATH_TYPE_KEYWORD
	#endif
#endif


//##############################################################################
//#                            RANLUX definitions                              #
//##############################################################################
#ifdef ENABLE_RANLUX
    #define RANLUXCL_LUX %(luxuaryFactor)s
    
    
    #ifdef ENABLE_DOUBLE
        #define RANLUXCL_SUPPORT_DOUBLE
        #define RAND_FLOAT_FUNCTION ranluxWrapper
        #define RANLUX_FUNCTION ranluxcl64
    #else
        #define RAND_FLOAT_FUNCTION ranluxWrapper
        #define RANLUX_FUNCTION ranluxcl32
    #endif
    
    
    #include "pyopencl-ranluxcl.cl"
    #define RAND_FUNCTION_ARG &ranluxclstate, &randCount, &random_temp
    #define RAND_FUNCTION_ARG_PTR ranluxclstate, randCount, random_temp
    #define RAND_INT_FUNCTION ranluxIntWrapper
    #define RAND_INT_FUNCTION_ARG &ranluxclstate, &randCount, &random_temp
#else
    #define RAND_FLOAT_FUNCTION randFloat
    #define RAND_FUNCTION_ARG &seed
    #define RAND_FUNCTION_ARG_PTR seedPtr
    #define RAND_INT_FUNCTION randInt
    #define RAND_INT_FUNCTION_ARG &seedG
#endif

//##############################################################################
//#                                 userCode                                   #
//##############################################################################
//Description: This is a placeholder for code that the user can write to help
//             specify operator and potential.
%(userCode)s


//##############################################################################
//#                                 OPERATOR                                   #
//##############################################################################
//Description: This returns the value of the operator for a specified point
//             in time (index in [0, N-1]).
#ifdef ENABLE_OPERATOR
inline void operator (DOF_ARGUMENT_DECL,FLOAT_TYPE opAccum[%(nbrOfOperators)s])
{
    %(operatorCode)s
}
#endif


//##############################################################################
//#                                 CORRELATOR                                 #
//##############################################################################
//Description: This returns the value of the operator for a specified point
//             in time (index in [0, N-1]).
#ifdef ENABLE_CORRELATOR
inline void correlator (DOF_ARGUMENT_DECL,FLOAT_TYPE corProd[%(nbrOfCorrelators)s])
{
    %(correlatorCode)s
}
#endif


//##############################################################################
//#                                 POTENTIAL                                  #
//##############################################################################
//Description: This returns the value of the potential for a specified
//             point in time (index in [0, N-1]).
inline FLOAT_TYPE potential(DOF_ARGUMENT_DECL)
{
    return %(potential)s;
}
//##############################################################################
//#                             ranlux_init_kernel                             #
//##############################################################################
//Description: Used to initialize and warmup the RANLUX PRNG. Run this as a
//             separate kernel instance before running the main kernel function.
#ifdef ENABLE_RANLUX
__kernel void ranlux_init_kernel(__global uint *seeds,
                                 __global ranluxcl_state_t *ranluxcltab)
{
    ranluxcl_initialization(*seeds, ranluxcltab);
}
//##############################################################################
//#                                   union                                    #
//##############################################################################
//Description: 
union ranlux_vector_union
{
    FLOAT_TYPE_VECTOR ranlux_vector;
    FLOAT_TYPE ranlux_array[4];
};
//##############################################################################
//#                               ranluxWrapper                                #
//##############################################################################
//Description: 
inline FLOAT_TYPE ranluxWrapper(ranluxcl_state_t *ranluxclstate,
                                uint *randCount,
                                union ranlux_vector_union *random_temp)
{
    if(((*randCount) & 3) == 0)
    {
        (*random_temp).ranlux_vector = RANLUX_FUNCTION (ranluxclstate);
    }
    
    return (*random_temp).ranlux_array[((*randCount)++) & 3];                                               // <---------
}
//##############################################################################
//#                             ranluxIntWrapper                               #
//##############################################################################
//Description: Used to produce four integers (uint4) with ranlux, ranging
//             between 0 and ranluxIntMax-1, using the same state buffer as the
//             other ranlux functions
inline uint ranluxIntWrapper(ranluxcl_state_t *ranluxclstate,
                            uint *randCount,
                            union ranlux_vector_union *random_temp)
{
    if(((*randCount) & 3) == 0)
    {
        (*random_temp).ranlux_vector = RANLUX_FUNCTION (ranluxclstate);
    }
    
    return (uint) ((*random_temp).ranlux_array[((*randCount)++) & 3] * ((FLOAT_TYPE) %(ranluxIntMax)s));     // <---------
}
#else // RANLUX NOT ENABLED
//##############################################################################
//#                                 XORSHIFT                                   #
//##############################################################################
//Description: This performs the necessary XORs to obtain a new random
//             uinteger w.
inline void xorshift (uint4 *seedPtr)
{
    uint t;
    t = (*seedPtr).x ^ ((*seedPtr).x << 11);
    (*seedPtr).xyz=(*seedPtr).yzw;
    (*seedPtr).w = ((*seedPtr).w ^ ((*seedPtr).w >> 19) ^ (t ^ (t >> 8)));
}
//##############################################################################
//#                               randFloat                                    #
//##############################################################################
//Description: This returns a random floating point number by dividing w with
//             UINT32_MAX (hardcoded).   
inline FLOAT_TYPE randFloat(uint4 *seedPtr)
{
    xorshift(seedPtr);
    return (*seedPtr).w * 2.328306437080797e-10;
}
//##############################################################################
//#                                randInt                                     #
//##############################################################################
//Description: This returns a random integer number through xorshift.
inline uint randInt(uint4 *seedPtr)
{
    xorshift(seedPtr);
    return (*seedPtr).w;
}
#endif // End of: RANLUX NOT ENABLED
//##############################################################################
//#                                kinEnergyEst                                #
//##############################################################################
//Description: This returns potential energy between the points.
inline FLOAT_TYPE kinEnergyEst (FLOAT_TYPE leftX, FLOAT_TYPE rightX)
{
    FLOAT_TYPE delta = leftX - rightX;
    return 0.5 * (delta * delta) * %(epsilon_inv2)s;
}
//##############################################################################
//#                                doBisectMove                                #
//##############################################################################
//Description: Does the bisection sampling algorithm for all DOF:s
//             and accepts/rejects the new path according to the
//             Metropolis algorithm.
#ifdef ENABLE_BISECTION

inline void doBisectMove (PATH_TYPE_KEYWORD FLOAT_TYPE *path,
#ifdef ENABLE_GLOBAL_OLD_PATH
                          __global FLOAT_TYPE *oldPath,
#endif
                         uint startPoint,
#ifdef ENABLE_RANLUX
                         ranluxcl_state_t *ranluxclstate,
                         uint *randCount,
                         union ranlux_vector_union *random_temp,
#else
                         uint4 *seedPtr,
#endif
                         uint *local_accepts, int *twoPow,
                         FLOAT_TYPE *sigmaN)
{
#ifndef ENABLE_GLOBAL_OLD_PATH
    // If path to be changed is store in local memory allocate
    // memory for it.
    FLOAT_TYPE oldPath[(%(2_POW_S)s -1) * %(DOF)s];
#endif

    // Save path that will later be changed. Also calculate action for the old
    // path.
    FLOAT_TYPE actionOld = 0.0;
    for (int i = 1; i < %(2_POW_S)s; i++)
    {
        for (int degree = 0; degree < %(DOF)s; degree++)
        {
            int Ndegree = degree * %(N)s;
            int p = i + Ndegree + startPoint;
            // Periodic boundary conditions.
            if (p >= Ndegree + %(N)s)
            {
                p -= %(N)s;
            }
            oldPath[(i - 1) * %(DOF)s + degree] = path[p];
        }

        int iCurr = i + startPoint;
        if (iCurr >= %(N)s)
        {
            iCurr -= %(N)s;
        }
        PATH_TYPE_KEYWORD FLOAT_TYPE *pathPointPtr = path + iCurr;
        actionOld += %(epsilon)s *potential(DOF_ARGUMENT_DATA);
    }

    // Loop over levels. Total of s levels. Starts from 1 to and including s.
    for (int n = 1; n <= %(S)s; n++)
    {
        // Loop over nodes for level n. Total of 2**(n-1) nodes to loop over.
        for (int i = 0; i < twoPow[n - 1]; i++)
        {
            // Move with normal dist with width sigmaN[n] = sqrt(epsilon * s**2 /
            // n**2)
            // First node to move will be located at 2**(s-n).
            // Distance between points are twice of that.
            // This gives:
            uint pointMoved =
                (1 + 2 * i) * twoPow[%(S)s - n] + startPoint;

            // The points to be averaged are pointMoved +- 2**(s-n) = i *
            // 2**(s-n+1) and (i+1) * 2**(s-n+1)
            uint pointLeft = i * twoPow[%(S)s - n + 1] + startPoint;
            uint pointRight = (i + 1) * twoPow[%(S)s - n + 1] + startPoint;

            // Periodic boundary conditions. None of the point values should be larger
            // than 2N.
            if (pointMoved >= %(N)s)
            {
                pointMoved = pointMoved - %(N)s;
            }
            if (pointLeft >= %(N)s)
            {
                pointLeft = pointLeft - %(N)s;
            }
            if (pointRight >= %(N)s)
            {
                pointRight = pointRight - %(N)s;
            }
            for (int degree = 0; degree < %(DOF)s; degree++)
            {
                int Ndegree = degree * %(N)s;
                // Need two random numbers for Box Muller method.
                FLOAT_TYPE u = RAND_FLOAT_FUNCTION(RAND_FUNCTION_ARG_PTR); //randFloat(seedPtr);
                FLOAT_TYPE v = RAND_FLOAT_FUNCTION(RAND_FUNCTION_ARG_PTR); //randFloat(seedPtr);
                FLOAT_TYPE rz;
                // Corner case if randFloat gives a zero.
                if (v == 0.0)
                {
                    rz = 0.0;
                }
                else
                {
                    rz = sigmaN[n] * sqrt (-2.0 * log (v))
                        * sin (2.0 * %(PI)s * u);
                }

                // Do the actual moving of node by averaging nodes to right and left
                // and use a normally distributed offset z.
                path[pointMoved + Ndegree] = 0.5 * (path[pointLeft + Ndegree] +
                                                     path[pointRight + Ndegree]) + rz;
            }
        }
    }

    // Calculate action for the new path
    FLOAT_TYPE actionNew = 0.0;
    for (int i = 1; i < %(2_POW_S)s; i++)
    {
        int iCurr = i + startPoint;
        if (iCurr >= %(N)s)
        {
            iCurr -= %(N)s;
        }
        PATH_TYPE_KEYWORD FLOAT_TYPE *pathPointPtr=path + iCurr;
        actionNew += %(epsilon)s *potential(DOF_ARGUMENT_DATA);
    }
    // Path rejected:
    if (exp(actionOld-actionNew) < RAND_FLOAT_FUNCTION(RAND_FUNCTION_ARG_PTR)) //randFloat(seedPtr))
    {
        // Revert to old path
        for (int i = 1; i < %(2_POW_S)s; i++)
        {
            for (int degree = 0; degree < %(DOF)s; degree++)
            {
                int Ndegree = degree * %(N)s;
                int p = i + Ndegree + startPoint;
                if (p >= Ndegree + %(N)s)
                {
                    p -= %(N)s;
                }
                path[p] = oldPath[(i - 1) * %(DOF)s + degree];
            }
        }
    }
    else
    {
        // Patch accepted, add one to accept counter.
        (*local_accepts)++;
    }
}
#endif
//##############################################################################
//#                                 shiftPathEnergyDiff                        #
//##############################################################################
//Description: Shifts all beads in an interval of a path by an offset and
//returns the difference in energy between the shifted and original path.
// TODO: This should really do the action checking itself and revert if needed.
#ifdef ENABLE_PATH_SHIFT
    inline FLOAT_TYPE
shiftPathEnergyDiff (PATH_TYPE_KEYWORD FLOAT_TYPE *x, uint degree, FLOAT_TYPE offset,
        uint left, uint right)
{
    if (left == right)
    {
        return 0.0;
    }
    FLOAT_TYPE kinDiff = 0.0;
    FLOAT_TYPE potDiff = 0.0;
    uint leftleft;
    uint rightright;
    PATH_TYPE_KEYWORD FLOAT_TYPE *xMinusNPart = x - %(N)s * degree;
    if (left < right)
    {
        rightright = right + 1;
        leftleft = left - 1;

        if (left == 0)
        {
            leftleft = %(N)s - 1;
        }
        if (right == %(N)s - 1)
        {
            rightright = 0;
        }

        kinDiff -= kinEnergyEst (x[left], x[leftleft]);
        kinDiff -= kinEnergyEst (x[right], x[rightright]);
        for (int i = left; i <= right; i++)
        {
	    PATH_TYPE_KEYWORD FLOAT_TYPE *pathPointPtr=xMinusNPart + i;
            potDiff -= potential (DOF_ARGUMENT_DATA);
            x[i] = x[i] + offset;
            potDiff += potential (DOF_ARGUMENT_DATA);
        }
        kinDiff += kinEnergyEst (x[left], x[leftleft]);
        kinDiff += kinEnergyEst (x[right], x[rightright]);

    }
    else
    {
        kinDiff -= kinEnergyEst (x[left], x[left - 1]);
        kinDiff -= kinEnergyEst (x[right], x[right + 1]);
        for (int i = left; i < %(N)s; i++)
        {
            PATH_TYPE_KEYWORD FLOAT_TYPE *pathPointPtr=xMinusNPart + i;
            potDiff -= potential (DOF_ARGUMENT_DATA);
            x[i] += offset;
            potDiff += potential (DOF_ARGUMENT_DATA);
        }

        for (int i = 0; i <= right; i++)
        {
            PATH_TYPE_KEYWORD FLOAT_TYPE *pathPointPtr=xMinusNPart + i;
            potDiff -= potential (DOF_ARGUMENT_DATA);
            x[i] += offset;
            potDiff += potential (DOF_ARGUMENT_DATA);
        }
        kinDiff += kinEnergyEst (x[left], x[left - 1]);
        kinDiff += kinEnergyEst (x[right], x[right + 1]);
    }
    return kinDiff + potDiff;

}
#endif
//##############################################################################
//#                                 Histogram                                  #
//##############################################################################
//Description: Function used to make a histogram representing the wavefunction
#ifdef ENABLE_BINS
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
    inline void
histogram (PATH_TYPE_KEYWORD FLOAT_TYPE *local_path, __global uint *binCounts)
{
    int ok = 1;
    int inc = 1;
    int binIndex = 0;

    for (int p = 0; p < %(DOF)s; p++)
    {
        if (local_path[%(N)s * p] < %(xMin)s ||
                local_path[%(N)s * p] >= %(xMax)s)
        {
            ok = 0;
            break;
        }

        binIndex +=
            (uint) ((local_path[%(N)s * p] - (%(xMin)s)) *
                    %(invBinSize)s)*inc;
        inc *= %(binsPerPart)s;
    }

    if (ok)
    {
        atom_inc(binCounts+binIndex);
    }
}
#endif

//##############################################################################
//#                                KERNEL CODE                                 #
//##############################################################################
//Description: The following code is the OpenCL kernel. This is the code that is
//             called from Python and takes all the input and delivers output.
__kernel void
metropolis (__global FLOAT_TYPE *paths
            ,__global uint *accepts
            ,__global uint *seeds
#ifdef ENABLE_GLOBAL_OLD_PATH
            ,__global FLOAT_TYPE *oldPaths
#endif
#ifdef ENABLE_OPERATOR
            ,__global FLOAT_TYPE *opMeans
#endif
#ifdef ENABLE_CORRELATOR
            ,__global FLOAT_TYPE *corMeans
#endif
#ifdef ENABLE_BINS
	        ,__global uint *binCounts
#endif
#ifdef ENABLE_RANLUX
            ,__global ranluxcl_state_t *ranluxcltab
#endif
        )
{
#ifdef ENABLE_PARALELLIZE_PATH
    uint walkerId = get_group_id(0) * %(nbrOfWalkersPerWorkGroup)s + get_local_id(1);       //<---------------
#endif
    uint threadId = get_global_id(0) + get_global_id(1) * get_global_size(0);
    
#ifdef ENABLE_RANLUX
    //Initialize the ranlux generator
    ranluxcl_initialization(*seeds, ranluxcltab);
    
    //ranluxclstate stores the state of the generator.
    ranluxcl_state_t ranluxclstate;

    //Download state into ranluxclstate struct.
    ranluxcl_download_seed(&ranluxclstate, ranluxcltab);
    
    union ranlux_vector_union random_temp;
    uint randCount = 0;
#else
    //This sets the seeds for the PRNG.
    uint4 seed,seedG;
    
    seed.x = seeds[threadId * 4 + 0];
    seed.y = seeds[threadId * 4 + 1];
    seed.z = seeds[threadId * 4 + 2];
    seed.w = seeds[threadId * 4 + 3];

    uint lastSeedPos = get_global_size(0) * get_global_size(1)*4;
    seedG.x = seeds[lastSeedPos + 0];
    seedG.y = seeds[lastSeedPos + 1];
    seedG.z = seeds[lastSeedPos + 2];
    seedG.w = seeds[lastSeedPos + 3];
#endif

    uint local_accepts=0;

#ifdef ENABLE_PATH_SHIFT
    uint mask = %(N)s - 1;
#endif

#ifdef ENABLE_GLOBAL_PATH
#ifdef ENABLE_PARALELLIZE_PATH
    PATH_TYPE_KEYWORD FLOAT_TYPE* local_path =paths + walkerId * %(pathSize)s;
#else
    PATH_TYPE_KEYWORD FLOAT_TYPE* local_path =paths + threadId * %(pathSize)s;
#endif
#else
    //This imports the path corresponding to this thread from the collection
    //of all paths stored in the field paths.
#ifdef ENABLE_PARALELLIZE_PATH
    PATH_TYPE_KEYWORD FLOAT_TYPE workGroupPath[%(pathSize)s * %(nbrOfWalkersPerWorkGroup)s];
    PATH_TYPE_KEYWORD FLOAT_TYPE *local_path = workGroupPath + get_local_id(1) * %(pathSize)s;

    for (uint i = 0; i < %(2_POW_S)s * %(DOF)s; i++)
        local_path[i + get_local_id(0) * %(2_POW_S)s * %(DOF)s] = paths[walkerId *
            %(pathSize)s + i + get_local_id(0) * %(2_POW_S)s * %(DOF)s];
#else
	PATH_TYPE_KEYWORD FLOAT_TYPE local_path[%(pathSize)s];
    for (uint i = 0; i < %(pathSize)s; i++)
        local_path[i] = paths[threadId * %(pathSize)s + i];
#endif
#endif

#ifdef ENABLE_BISECTION
    // Create the stdev array for bisection algorithm and
    // an array with powers of two.
    int twoPow[%(S)s + 1];
    FLOAT_TYPE sigmaN[%(S)s + 1];
    int Ns = %(2_POW_S)s;
    twoPow[0] = 1;

    sigmaN[0] = ((FLOAT_TYPE) Ns) * %(epsilon)s;
    for (int i = 1; i <= %(S)s; i++)
    {
        twoPow[i] = twoPow[i - 1] * 2;
        sigmaN[i] = sqrt (((FLOAT_TYPE) Ns * %(epsilon)s / (2.0 * (FLOAT_TYPE) twoPow[i])));
    }
#endif

#ifdef ENABLE_OPERATOR
    FLOAT_TYPE opAccum[%(nbrOfOperators)s]={%(nbrOfOperatorsZeros)s};
#endif

    // The following loops are responsible for creating Metropolis
    // samples and measuring the operator.

    for (uint i = 0; i < %(operatorRuns)s; i++)
    {

#ifdef ENABLE_PARALELLIZE_PATH
#ifdef ENABLE_GLOBAL_PATH
        barrier(CLK_GLOBAL_MEM_FENCE);
#else
        barrier(CLK_LOCAL_MEM_FENCE);
#endif
#endif
        for(uint j = 0; j < %(metroStepsPerOperatorRun)s; j++)
        {
#ifdef ENABLE_PATH_SHIFT
            uint degree = (%(metroStepsPerOperatorRun)s*i+j) %% %(DOF)s;
            // Choose the interval of the path to shift
            //xorshift (&seedG);
            //uint left =seedG.w & mask;
            //xorshift (&seedG);
            //uint right = seedG.w & mask;
            uint left = RAND_INT_FUNCTION(RAND_INT_FUNCTION_ARG) & mask;                                    //<-------------
            uint right = RAND_INT_FUNCTION(RAND_INT_FUNCTION_ARG) & mask;
            FLOAT_TYPE offset =
                %(PSAlpha)s * (2.0 * RAND_FLOAT_FUNCTION(RAND_FUNCTION_ARG) - 1.0); //randFloat (&seed) - 1.0);

            // Revert path shift if rejected by Metropolis step
            if (exp (-%(epsilon)s *  shiftPathEnergyDiff (local_path + (%(N)s * degree),
                            degree, offset, left, right)) < RAND_FLOAT_FUNCTION(RAND_FUNCTION_ARG)) //randFloat (&seed))
            {
                shiftPathEnergyDiff (local_path + (%(N)s * degree),
                        degree, -offset, left, right);
            }
            else
            {
                local_accepts++;
            }
#endif

#ifdef ENABLE_SINGLE_NODE_MOVE
// TODO: Make into its own function
            uint modPoint = (i * %(metroStepsPerOperatorRun)s + j)%% %(N)s;
            uint modPathPoint = (i * %(metroStepsPerOperatorRun)s+j)%% (%(N)s * %(DOF)s);
            //Arrange with the periodic boundary conditions.
            uint rightIndex = modPathPoint + 1;
            uint leftIndex = modPathPoint - 1;
            if (modPoint == 0)
            {
                leftIndex += %(N)s;
            }
            else if (modPoint == (%(N)s - 1))
            {
                rightIndex -= %(N)s;
            }

            FLOAT_TYPE oldX = local_path[modPathPoint];
            FLOAT_TYPE modX = oldX + (2.0 *RAND_FLOAT_FUNCTION(RAND_FUNCTION_ARG)-1.0)*%(alpha)s;
            //FLOAT_TYPE modX = oldX + (2.0 *randFloat(&seed)-1.0)*%(alpha)s;

            //Calculate the difference in energy (action) for the new path
            //compared to the old, stored in diffE.
            FLOAT_TYPE diffE = kinEnergyEst(modX, local_path[rightIndex])
                + kinEnergyEst(modX, local_path[leftIndex])
                - kinEnergyEst(local_path[modPathPoint], local_path[leftIndex])
                - kinEnergyEst(local_path[modPathPoint],
                        local_path[rightIndex]);

            PATH_TYPE_KEYWORD FLOAT_TYPE *pathPointPtr=local_path + modPoint;
            diffE -= potential(DOF_ARGUMENT_DATA);
            local_path[modPathPoint] = modX;
            diffE += potential(DOF_ARGUMENT_DATA);
            local_path[modPathPoint] = oldX;

            //Determine whether or not to accept the change in the path.
            if (native_exp(-%(epsilon)s*diffE) > RAND_FLOAT_FUNCTION(RAND_FUNCTION_ARG)) //randFloat(&seed))
            {
                local_path[modPathPoint] = modX;
                local_accepts++;
            }
#endif

#ifdef ENABLE_BISECTION
            //xorshift(&seedG);
            doBisectMove(local_path,
#ifdef ENABLE_GLOBAL_OLD_PATH
            oldPaths+threadId*(%(2_POW_S)s-1)*%(DOF)s,
#endif
#ifdef ENABLE_PARALELLIZE_PATH
             
             ((RAND_INT_FUNCTION(RAND_INT_FUNCTION_ARG) & ((uint)(%(2_POW_S)s - 1))) + %(2_POW_S)s*get_local_id(0))%% %(N)s,
             //((seedG.w & ((uint)(%(2_POW_S)s - 1))) + %(2_POW_S)s*get_local_id(0))%% %(N)s,
#else
             (RAND_INT_FUNCTION(RAND_INT_FUNCTION_ARG) & ((uint)(%(N)s - 1)))%% %(N)s,
             //(seedG.w & ((uint)(%(N)s - 1)))%% %(N)s,
#endif
             RAND_FUNCTION_ARG,&local_accepts, twoPow, sigmaN);
#endif

#ifdef ENABLE_PARALELLIZE_PATH
#ifdef ENABLE_GLOBAL_PATH
            barrier(CLK_GLOBAL_MEM_FENCE);
#else
            barrier(CLK_LOCAL_MEM_FENCE);
#endif
#endif
        }//end of metroStepsPerOperatorRun loop


#if (defined ENABLE_OPERATOR) || (defined ENABLE_BINS) || (defined ENABLE_CORRELATOR)
#ifdef ENABLE_PARALELLIZE_PATH
        for(uint timeSlice = get_local_id(0) * %(2_POW_S)s;
                timeSlice<(get_local_id(0)+1) * %(2_POW_S)s; timeSlice++)
#else
        for(uint timeSlice = 0; timeSlice < %(N)s; timeSlice++)
#endif
        {
#ifdef ENABLE_OPERATOR
	    PATH_TYPE_KEYWORD FLOAT_TYPE *pathPointPtr = local_path + timeSlice;
	    operator(DOF_ARGUMENT_DATA, opAccum);
#endif
#ifdef ENABLE_BINS
            histogram (local_path + timeSlice,binCounts);
#endif
#ifdef ENABLE_CORRELATOR
#ifdef ENABLE_PARALELLIZE_PATH
	    for(uint corTI = 0; corTI< %(N)s / 2; corTI++)
            {
	        uint corT = corTI + get_local_id(0);
		if(corT >= %(N)s/2)
                {
		    corT-=%(N)s/2;
                }
            
#else
	    for(uint corT=0;corT < %(N)s / 2; corT++)
            {
#endif
	        uint timeSlice2;

		if(timeSlice+corT < %(N)s)
                {
                    timeSlice2 = timeSlice + corT;
                }
	        else
                {
		    timeSlice2 = timeSlice + corT - %(N)s;
                }

    	            FLOAT_TYPE corProd[%(nbrOfCorrelators)s] = {%(nbrOfCorrelatorsOnes)s};
                    PATH_TYPE_KEYWORD FLOAT_TYPE *pathPointPtr = local_path + timeSlice;
                    correlator(DOF_ARGUMENT_DATA, corProd);

                    pathPointPtr = local_path + timeSlice2;

                    correlator(DOF_ARGUMENT_DATA, corProd);
                    for(uint corI = 0; corI < %(nbrOfCorrelators)s; corI++)
                    {
#ifdef ENABLE_PARALELLIZE_PATH
                        corMeans[walkerId * %(nbrOfCorrelators)s * %(N)s/2 +
                                 corI * %(N)s/2 + corT] += corProd[corI];
#else
                        corMeans[threadId * %(nbrOfCorrelators)s * %(N)s/2 +
                                 corI * %(N)s/2 + corT] += corProd[corI];
#endif
		    }
#ifdef ENABLE_PARALELLIZE_PATH
		    barrier(CLK_GLOBAL_MEM_FENCE);
#endif
         	}
#endif
	    }
#endif

	}//end of operatorRuns loop

#ifdef ENABLE_OPERATOR
	for(int op = 0; op < %(nbrOfOperators)s; op++)
    	opMeans[op + threadId * %(nbrOfOperators)s] = opAccum[op] * %(opNorm)s;
#endif

#ifndef ENABLE_GLOBAL_PATH
    //Store the last path back in the paths variable for use in the
    //next kernel run.
#ifdef ENABLE_PARALELLIZE_PATH
    for (uint i = 0; i < %(2_POW_S)s * %(DOF)s; i++)
    {
        paths[walkerId * %(pathSize)s + i + get_local_id(0) * %(2_POW_S)s * %(DOF)s]
              = local_path[i + get_local_id(0) * %(2_POW_S)s * %(DOF)s];
    }
#else
    for (uint i = 0; i < %(pathSize)s; i++)
    {
        paths[threadId * %(pathSize)s + i] = local_path[i];
    }
#endif
#endif

    accepts[threadId] = local_accepts;

#ifdef ENABLE_RANLUX
    //Upload ranlux state for the consequitve runs.
    ranluxcl_upload_seed(&ranluxclstate, ranluxcltab);
#else
    //Store the current state of the Xorshift PRNG for use in the next
    //kernel run.
    seeds[threadId * 4] = seed.x;
    seeds[threadId * 4 + 1] = seed.y;
    seeds[threadId * 4 + 2] = seed.z;
    seeds[threadId * 4 + 3] = seed.w;
#endif
}

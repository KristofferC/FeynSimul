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
inline void operator (DOF_ARGUMENT_DECL,float opAccum[%(nbrOfOperators)s])
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
inline void correlator (DOF_ARGUMENT_DECL,float corProd[%(nbrOfCorrelators)s])
{
    %(correlatorCode)s
}
#endif


//##############################################################################
//#                                 POTENTIAL                                  #
//##############################################################################
//Description: This returns the value of the potential for a specified
//             point in time (index in [0, N-1]).
inline float potential(DOF_ARGUMENT_DECL)
{
    return %(potential)s;
}

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
inline float
randFloat(uint4 *seedPtr)
{
    xorshift(seedPtr);
    return (*seedPtr).w * 2.328306437080797e-10;
}


//##############################################################################
//#                                kinEnergyEst                                #
//##############################################################################
//Description: This returns potential energy between the points.
inline float kinEnergyEst (float leftX, float rightX)
{
    float delta = leftX - rightX;
    return 0.5f * (delta * delta) * %(epsilon_inv2)s;
}

//##############################################################################
//#                                doBisectMove                                #
//##############################################################################
//Description: Does the bisection sampling algorithm for all DOF:s
//             and accepts/rejects the new path according to the
//             Metropolis algorithm.
#ifdef ENABLE_BISECTION

inline void doBisectMove (PATH_TYPE_KEYWORD float *path,
#ifdef ENABLE_GLOBAL_OLD_PATH
                          __global float *oldPath,
#endif
                         uint startPoint, uint4 *seedPtr,
                         uint *local_accepts, int *twoPow,
                         float *sigmaN)
{
#ifndef ENABLE_GLOBAL_OLD_PATH
    // If path to be changed is store in local memory allocate
    // memory for it.
    float oldPath[(%(2_POW_S)s -1) * %(DOF)s];
#endif

    // Save path that will later be changed. Also calculate action for the old
    // path.
    float actionOld = 0.0f;
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
        PATH_TYPE_KEYWORD float *pathPointPtr = path + iCurr;
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
                float u = randFloat(seedPtr);
                float v = randFloat(seedPtr);
                float rz;
                // Corner case if randFloat gives a zero.
                if (v == 0.0f)
                {
                    rz = 0.0f;
                }
                else
                {
                    rz = sigmaN[n] * native_sqrt (-2.0f * log (v))
                        * native_sin (2.0f * %(PI)s * u);
                }

                // Do the actual moving of node by averaging nodes to right and left
                // and use a normally distributed offset z.
                path[pointMoved + Ndegree] = 0.5f * (path[pointLeft + Ndegree] +
                                                     path[pointRight + Ndegree]) + rz;
            }
        }
    }

    // Calculate action for the new path
    float actionNew = 0.0f;
    for (int i = 1; i < %(2_POW_S)s; i++)
    {
        int iCurr = i + startPoint;
        if (iCurr >= %(N)s)
        {
            iCurr -= %(N)s;
        }
        PATH_TYPE_KEYWORD float *pathPointPtr=path + iCurr;
        actionNew += %(epsilon)s *potential(DOF_ARGUMENT_DATA);
    }
    // Path rejected:
    if (exp(actionOld-actionNew) < randFloat(seedPtr))
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
    inline float
shiftPathEnergyDiff (PATH_TYPE_KEYWORD float *x, uint degree, float offset,
        uint left, uint right)
{
    if (left == right)
    {
        return 0.0f;
    }
    float kinDiff = 0.0f;
    float potDiff = 0.0f;
    uint leftleft;
    uint rightright;
    PATH_TYPE_KEYWORD float *xMinusNPart = x - %(N)s * degree;
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
	    PATH_TYPE_KEYWORD float *pathPointPtr=xMinusNPart + i;
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
            PATH_TYPE_KEYWORD float *pathPointPtr=xMinusNPart + i;
            potDiff -= potential (DOF_ARGUMENT_DATA);
            x[i] += offset;
            potDiff += potential (DOF_ARGUMENT_DATA);
        }

        for (int i = 0; i <= right; i++)
        {
            PATH_TYPE_KEYWORD float *pathPointPtr=xMinusNPart + i;
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
histogram (PATH_TYPE_KEYWORD float *local_path, __global uint *binCounts)
{
    int ok = 1;
    int inc = 1;
    int binIndex = 0;

    for (int p = 0; p < %(DOF)s; p++)
    {
        if (local_path[%(N)s * p] < %(xmin)s ||
                local_path[%(N)s * p] >= %(xmax)s)
        {
            ok = 0;
            break;
        }

        binIndex +=
            (uint) ((local_path[%(N)s * p] - (%(xmin)s)) *
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
metropolis (__global float *paths
            ,__global uint *accepts
            ,__global uint *seeds
#ifdef ENABLE_GLOBAL_OLD_PATH
            ,__global float *oldPaths
#endif
#ifdef ENABLE_OPERATOR
           ,__global float *opMeans
#endif
#ifdef ENABLE_CORRELATOR
            ,__global float *corMeans
#endif
#ifdef ENABLE_BINS
	    ,__global uint *binCounts
#endif
        )
{
#ifdef ENABLE_PARALELLIZE_PATH
    uint walkerId = get_group_id(0) * %(nbrOfWalkersPerWorkGroup)s + get_local_id(1);
#endif
    uint threadId = get_global_id(0) + get_global_id(1) * get_global_size(0);
    //This sets the seeds for the Xorshift PRNG.
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

    uint local_accepts=0;

#ifdef ENABLE_PATH_SHIFT
    uint mask = %(N)s - 1;
#endif

#ifdef ENABLE_GLOBAL_PATH
#ifdef ENABLE_PARALELLIZE_PATH
    PATH_TYPE_KEYWORD float* local_path =paths + walkerId * %(pathSize)s;
#else
    PATH_TYPE_KEYWORD float* local_path =paths + threadId * %(pathSize)s;
#endif
#else
    //This imports the path corresponding to this thread from the collection
    //of all paths stored in the field paths.
#ifdef ENABLE_PARALELLIZE_PATH
    PATH_TYPE_KEYWORD float workGroupPath[%(pathSize)s * %(nbrOfWalkersPerWorkGroup)s];
    PATH_TYPE_KEYWORD float *local_path = workGroupPath + get_local_id(1) * %(pathSize)s;

    for (uint i = 0; i < %(2_POW_S)s * %(DOF)s; i++)
        local_path[i + get_local_id(0) * %(2_POW_S)s * %(DOF)s] = paths[walkerId *
            %(pathSize)s + i + get_local_id(0) * %(2_POW_S)s * %(DOF)s];
#else
	PATH_TYPE_KEYWORD float local_path[%(pathSize)s];
    for (uint i = 0; i < %(pathSize)s; i++)
        local_path[i] = paths[threadId * %(pathSize)s + i];
#endif
#endif

#ifdef ENABLE_BISECTION
    // Create the stdev array for bisection algorithm and
    // an array with powers of two.
    int twoPow[%(S)s + 1];
    float sigmaN[%(S)s + 1];
    int Ns = %(2_POW_S)s;
    twoPow[0] = 1;

    sigmaN[0] = ((float) Ns) * %(epsilon)s;
    for (int i = 1; i <= %(S)s; i++)
    {
        twoPow[i] = twoPow[i - 1] * 2;
        sigmaN[i] = sqrt (((float) Ns * %(epsilon)s / (2.0f * (float) twoPow[i])));
    }
#endif

#ifdef ENABLE_OPERATOR
    float opAccum[%(nbrOfOperators)s]={%(nbrOfOperatorsZeros)s};
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
            xorshift (&seedG);
            uint left =seedG.w & mask;
            xorshift (&seedG);
            uint right = seedG.w & mask;
            float offset =
                %(PSAlpha)s * (2.0f * randFloat (&seed) - 1.0f);

            // Revert path shift if rejected by Metropolis step
            if (exp (-%(epsilon)s *  shiftPathEnergyDiff (local_path + (%(N)s * degree),
                            degree, offset, left, right)) < randFloat (&seed))
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

            float oldX = local_path[modPathPoint];
            float modX = oldX + (2.0f *randFloat(&seed)-1.0f)*%(alpha)s;

            //Calculate the difference in energy (action) for the new path
            //compared to the old, stored in diffE.
            float diffE = kinEnergyEst(modX, local_path[rightIndex])
                + kinEnergyEst(modX, local_path[leftIndex])
                - kinEnergyEst(local_path[modPathPoint], local_path[leftIndex])
                - kinEnergyEst(local_path[modPathPoint],
                        local_path[rightIndex]);

            PATH_TYPE_KEYWORD float *pathPointPtr=local_path + modPoint;
            diffE -= potential(DOF_ARGUMENT_DATA);
            local_path[modPathPoint] = modX;
            diffE += potential(DOF_ARGUMENT_DATA);
            local_path[modPathPoint] = oldX;

            //Determine whether or not to accept the change in the path.
            if (native_exp(-%(epsilon)s*diffE) > randFloat(&seed))
            {
                local_path[modPathPoint] = modX;
                local_accepts++;
            }
#endif

#ifdef ENABLE_BISECTION
            xorshift(&seedG);
            doBisectMove(local_path,
#ifdef ENABLE_GLOBAL_OLD_PATH
            oldPaths+threadId*(%(2_POW_S)s-1)*%(DOF)s,
#endif
#ifdef ENABLE_PARALELLIZE_PATH
             ((seedG.w & ((uint)(%(2_POW_S)s - 1))) + %(2_POW_S)s*get_local_id(0))%% %(N)s,
#else
             (seedG.w & ((uint)(%(N)s - 1)))%% %(N)s,
#endif
             &seed,&local_accepts, twoPow, sigmaN);
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
	    PATH_TYPE_KEYWORD float *pathPointPtr = local_path + timeSlice;
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

    	            float corProd[%(nbrOfCorrelators)s] = {%(nbrOfCorrelatorsOnes)s};
                    PATH_TYPE_KEYWORD float *pathPointPtr = local_path + timeSlice;
                    correlator(DOF_ARGUMENT_DATA, corProd);

                    pathPointPtr = local_path + timeSlice2;

                    correlator(DOF_ARGUMENT_DATA, corProd);
                    for(uint corI = 0; corI < %(nbrOfCorrelators)s; corI++)
                    {
#ifdef ENABLE_PARALELLIZE_PATH
                        corMeans[walkerId * %(nbrOfCorrelators)s * %(N)s/2 +
                                 corI * %(N)s/2  +corT]+=corProd[corI];
#else
                        corMeans[threadId*%(nbrOfCorrelators)s*%(N)s/2  +
                                 corI*%(N)s/2  +corT]+=corProd[corI];
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
              = local_path[i + get_local_id(0)*%(2_POW_S)s* %(DOF)s];
    }
#else
    for (uint i = 0; i < %(pathSize)s; i++)
    {
        paths[threadId * %(pathSize)s + i] = local_path[i];
    }
#endif
#endif

    accepts[threadId] = local_accepts;

    //Store the current state of the Xorshift PRNG for use in the next
    //kernel run.
    seeds[threadId * 4] = seed.x;
    seeds[threadId * 4 + 1] = seed.y;
    seeds[threadId * 4 + 2] = seed.z;
    seeds[threadId * 4 + 3] = seed.w;
}

#############################################################################
# ___________________________________________________________________________
# Authors:     Kristoffer Carlsson - kricarl@student.chalmers.se
#              Patric Holmvall - holmvall@student.chalmers.se
#              Petter Saterskog - spetter@student.chalmers.se
#
# Date:        2011-07-29
#

"""

"""
__docformat__ = 'epytext en'


import sys
import numpy as np
from time import time
from copy import copy
from collections import defaultdict
import os

import pyopencl as cl
import pyopencl.array
import pyopencl.clrandom
import pyopencl.clmath

from pimc_kernel import PIMCKernel


def loadKernel(system, RP):
    """
    Loads the kernel to the GPU and readies it to be run.

    The kernel code from kernel.src is read and modified with the necessary
    defines and dictionary replacements. The binary program that will be
    sent to the GPU is built and memory allocations on the GPU are done.

    @type system: L{physical_systems} class
    @param system: The physical system that is considered.
    @type RP: L{run_params.RunParams}
    @param RP: The run parameters that have been set for the run.
    @rtype: L{pimc_kernel.PIMCKernel}
    @return: Class containing the stuff needed by the function L{runKernel}.
    """

    pimcKernel = PIMCKernel()

    if not (RP.nbrOfWalkers / RP.nbrOfWalkersPerWorkGroup *
           RP.nbrOfWalkersPerWorkGroup == RP.nbrOfWalkers):
        raise Exception("The number of walkers must be a multiple "
                        "of the number of walkers per work group.")

    if  RP.nbrOfWalkers < RP.nbrOfWalkersPerWorkGroup:
        raise Exception("The number of walkers per group must be less "
                        "than or equal to the total number of walkers")

    if  not (RP.N != 0 and ((RP.N & (RP.N - 1)) == 0)):
        raise Exception("Path node number, N, must be a power of 2.")

    if RP.enableBisection:
        if 2 ** RP.S > RP.N:
            raise Exception("2^S must be less than or equal to N.")

    if RP.enableParallelizePath:
        if RP.enableSingleNodeMove:
            raise Exception("Single node move not possible when "
                            "enableParallelizePath is True.")
        if RP.enablePathShift:
            raise Exception("Shift path not possible when "
                            "enableParallelizePath is True.")

        nbrOfThreadsPerWalker = RP.N / (2 ** RP.S)
        nbrOfWorkgroups = RP.nbrOfWalkers / RP.nbrOfWalkersPerWorkGroup

        pimcKernel.localSize = (nbrOfThreadsPerWalker,
                                RP.nbrOfWalkersPerWorkGroup)
        pimcKernel.globalSize = (nbrOfWorkgroups * nbrOfThreadsPerWalker,
                                 RP.nbrOfWalkersPerWorkGroup)
        pimcKernel.nbrOfThreads = RP.nbrOfWalkers * nbrOfThreadsPerWalker
    else:
        pimcKernel.localSize = None
        pimcKernel.globalSize = (RP.nbrOfWalkers,)
        pimcKernel.nbrOfThreads = RP.nbrOfWalkers

    defines = ""

    if RP.enableBins:
        defines += "#define ENABLE_BINS\n"

    if RP.enableOperator:
        defines += "#define ENABLE_OPERATOR\n"

    if RP.enableCorrelator:
        defines += "#define ENABLE_CORRELATOR\n"

    if RP.enablePathShift:
        defines += "#define ENABLE_PATH_SHIFT\n"

    if RP.enableBisection:
        defines += "#define ENABLE_BISECTION\n"

    if RP.enableSingleNodeMove:
        defines += "#define ENABLE_SINGLE_NODE_MOVE\n"

    if RP.enableGlobalPath:
        defines += "#define ENABLE_GLOBAL_PATH\n"

    if RP.enableGlobalOldPath:
        if not RP.enableBisection:
            raise Exception("Old Path is not used when "
                            "bisection is not enabled")
        defines += "#define ENABLE_GLOBAL_OLD_PATH\n"

    if RP.enableParallelizePath:
        defines += "#define ENABLE_PARALELLIZE_PATH\n"

    defines += "#define DOF_ARGUMENT_DECL "
    for i in range(system.DOF):
        defines += "float x" + str(i + 1)
        if (i != system.DOF - 1):
            defines += ","
    defines += "\n"

    defines += "#define DOF_ARGUMENT_DATA "
    for i in range(system.DOF):
        defines += "pathPointPtr[" + str(i) + " * " + str(RP.N) + "]"
        if (i != system.DOF - 1):
            defines += ","
    defines += "\n"

    if RP.enableOperator:
        operatorCode = ""
        nbrOfOperatorsZeros = ""
        for i in range(0, len(RP.operators)):
            operatorCode += ("opAccum[" + str(i) + "] += "
                             + RP.operators[i] + ";\n")
            if i != 0:
                nbrOfOperatorsZeros += ","
            nbrOfOperatorsZeros += "0.0f"

    if RP.enableCorrelator:
        correlatorCode = ""
        nbrOfCorrelatorsOnes = ""
        for i in range(0, len(RP.correlators)):
            correlatorCode += ("corProd[" + str(i) + "] *= "
                               + RP.correlators[i] + ";\n")
            if i != 0:
                nbrOfCorrelatorsOnes += ","
            nbrOfCorrelatorsOnes += "1.0f"

    class DictWithDefault(defaultdict):
        def __missing__(self, key):
            return key + str(" is not defined")

    # Replacements that will be done in the kernel code file will be stored in
    # this dictionary.
    replacements = DictWithDefault()

    replacements['nbrOfWalkersPerWorkGroup'] = RP.nbrOfWalkersPerWorkGroup
    replacements['N'] = '%d' % RP.N
    replacements['potential'] = system.potential
    replacements['epsilon'] = '%ef' % (RP.beta / float(RP.N))
    replacements['epsilon_inv2'] = '%ef' % ((float(RP.N) / RP.beta) ** 2)
    replacements['pathSize'] = '%d' % (system.DOF * RP.N)
    replacements['defines'] = defines
    replacements['userCode'] = system.userCode
    replacements['DOF'] = '%d' % system.DOF
    replacements['operatorRuns'] = '%d' % RP.operatorRuns
    replacements['metroStepsPerOperatorRun'] = ('%d'
            % RP.metroStepsPerOperatorRun)

    if RP.enableOperator:
        replacements['operatorCode'] = operatorCode
        replacements['nbrOfOperators'] = '%d' % len(RP.operators)
        replacements['nbrOfOperatorsZeros'] = nbrOfOperatorsZeros
        replacements['opNorm'] = '%ef' % (1.0 / float(RP.operatorRuns * RP.N
                                * RP.nbrOfWalkers / pimcKernel.nbrOfThreads))

    if RP.enableCorrelator:
        replacements['correlatorCode'] = correlatorCode
        replacements['nbrOfCorrelators'] = '%d' % len(RP.correlators)
        replacements['nbrOfCorrelatorsOnes'] = nbrOfCorrelatorsOnes

    if RP.enableBisection:
        replacements['PI'] = '%ef' % 3.141592653589793
        replacements['2_POW_S'] = '%d' % 2 ** RP.S
        replacements['S'] = '%d' % RP.S

    if RP.enableSingleNodeMove:
        replacements['alpha'] = '%ef' % RP.alpha

    if RP.enablePathShift:
        replacements['PSAlpha'] = '%ef' % RP.PSAlpha

    if RP.enableBins:
        replacements['xmin'] = '%ef' % RP.xmin
        replacements['xmax'] = '%ef' % RP.xmax
        replacements['binsPerPart'] = '%d' % RP.binResolutionPerDOF
        replacements['nbrOfBins'] = '%d' % RP.binResolutionPerDOF ** system.DOF
        replacements['invBinSize'] = '%ef' % (float(RP.binResolutionPerDOF) /
                                     float(RP.xmax - RP.xmin))

    #Import kernel code and paste 'replacements' into it
    kernelCode_r = open(os.path.dirname(__file__) +
            '/GPUsrc/kernel.c', 'r').read()
    kernelCode = kernelCode_r % replacements
    """ String containing the kernel code that will be sent to the GPU. """

    printCleaned = False
    if printCleaned:
        import commands
        preprocessedCode = commands.getstatusoutput('echo "' +
                kernelCode + '" | cpp')[1]
        cleanedPreprocessedCode = ""
        for i in preprocessedCode.splitlines():
            if len(i) > 0:
                if i[0] != '#':
                    cleanedPreprocessedCode += i + '\n'
        print cleanedPreprocessedCode

    #Create the OpenCL context and command queue
    pimcKernel.ctx = cl.create_some_context()
    queueProperties = cl.command_queue_properties.PROFILING_ENABLE
    pimcKernel.queue = cl.CommandQueue(pimcKernel.ctx,
                                       properties=queueProperties)

    #Build the program and identify metropolis as the kernel
    pimcKernel.prg = (cl.Program(pimcKernel.ctx, kernelCode)
                     .build(options="-cl-fast-relaxed-math"))
    pimcKernel.kernel = pimcKernel.prg.metropolis

    #Initial paths are created (the initial path vector is filled with zeros,
    #meaning no movement of the particles)
    try:
        pimcKernel.paths = cl.array.zeros(pimcKernel.queue,
                          (RP.nbrOfWalkers, RP.N * system.DOF),
                          np.float32)

        #Buffer for storing number of accepted values and
        #seeds for the xorshfitPRNG
        pimcKernel.accepts = cl.array.zeros(pimcKernel.queue,
                (pimcKernel.nbrOfThreads, ), np.uint32)

        #np.random.seed(0)
        pimcKernel.seeds = cl.array.to_device(pimcKernel.queue,
                         (np.random.randint(0, high = 2 ** 31 -
                          size = (pimcKernel.nbrOfThreads + 1, 4))
                          ).astype(np.uint32))
        if RP.enableOperator:
            #pyopencl.array objects are created for storing
            #the calculated operator means from each thread
            pimcKernel.operatorValues = cl.array.zeros(pimcKernel.queue,
                    pimcKernel.nbrOfThreads * len(RP.operators), np.float32)

        if RP.enableCorrelator:
            #pyopencl.array objects are created for storing
            #the calculated operator means from each thread
            pimcKernel.correlatorValues = cl.array.zeros(pimcKernel.queue,
                    (RP.nbrOfWalkers, len(RP.correlators), RP.N / 2),
                    np.float32)

        if RP.enableGlobalOldPath:
            pimcKernel.oldPath = cl.array.zeros(pimcKernel.queue,
                                (pimcKernel.nbrOfThreads,
                                (2 ** RP.S - 1) * system.DOF), np.float32)

        #A buffer for the multidimensional array for storing
        #the resulting number of bin counts
        if RP.enableBins:
            binTouple = (RP.binResolutionPerDOF,)
            for i in range(1, system.DOF):
                binTouple += (RP.binResolutionPerDOF,)
            pimcKernel.binCounts = cl.array.zeros(pimcKernel.queue,
                    binTouple, np.uint32)

    except pyopencl.MemoryError:
        raise Exception("Unable to allocate global "
                        "memory on device, out of memory?")

    pimcKernel.runParams = copy(RP)
    pimcKernel.system = copy(system)
    pimcKernel.nbrOfOperators = len(RP.operators)

    #Return the kernel object, context, command queue and
    #pyopencl.array objects
    return pimcKernel


class RunKernelResults:
    """
    RunKernelResults is a class where results from L{runKernel} are dumped.

    After a kernel run the following variables are available:

    Always:
      - B{runTime} (float) E{-} The time in seconds it took to execute
        the kernel.
      - B{acceptanceRate} (float)- E{-} Average of how often path
        changes where accepted or rejected in the Metropolis step.
    If B{getOperator} in L{run_params.RunParams} is enabled:
      - B{operatorMean} (float) E{-} The average value of the operators
        calculated by the threads.
      - B{operatorStandardError} (float) E{-} The standard error in
        the operatorMean.
    If B{binsEnabled} and B{returnBinCounts} in L{run_params.RunParams} are
    enabled:
      - B{binCountNormed} (numpy array) E{-} Calculated probability density
        in the interval specified in the run parameters.
    If B{returnPaths} in L{run_params} is enabled:
      - B{paths} (numpy array) E{-} All the current active
        paths for the threads.
    """
    pass


def runKernel(pimcKernel):
    """
    Runs the kernel with the system and the run parameters.

    This will start the kernel on the GPU. When the GPU is done desired
    data is fetched.

    @type pimcKernel: L{pimc_kernel.PIMCKernel}
    @param pimcKernel: The class that is returned from L{loadKernel}.
    @rtype: L{RunKernelResults}
    @return: Contains the results of running the kernel.
    """
    runParams = pimcKernel.runParams

    RKR = RunKernelResults()

    #Run kernel. The buffers for all pyopencl.array objects are passed
    #by the .data parameter.

    args = [pimcKernel.paths.data,
            pimcKernel.accepts.data,
            pimcKernel.seeds.data]

    if runParams.enableGlobalOldPath:
        args += [pimcKernel.oldPath.data]

    if runParams.enableOperator:
        args += [pimcKernel.operatorValues.data]

    if runParams.enableCorrelator:
        args += [pimcKernel.correlatorValues.data]

    if runParams.enableBins:
        args += [pimcKernel.binCounts.data]

    #print(pimcKernel.globalSize)#(nbrOfWorkgroups* nbrOfThreadsPerWalker, RP.nbrOfWalkersPerWorkGroup)
    #print(pimcKernel.localSize)#(nbrOfThreadsPerWalker,RP.nbrOfWalkersPerWorkGroup)
    #print(pimcKernel.getStats())

    startTime = time()
    kernelObj = pimcKernel.kernel(pimcKernel.queue,
                                      pimcKernel.globalSize,
                                      pimcKernel.localSize,
                                      *args)

    #Wait until the threads have finished and then calculate total run time
    try:
        kernelObj.wait()
    except pyopencl.RuntimeError:
        if time() - startTime > 5.0:
            raise Exception("Kernel runtime error. 5 seconds passed "
                            "and kernel aborted.")
        else:
            raise

    RKR.runTime = 1e-9 * (kernelObj.profile.end - kernelObj.profile.start)

    #Fetch the operatorValues and acceptanceRate pyopencl.array objects from
    #the graphics card to the RAM.
    RKR.acceptanceRate = (pimcKernel.accepts.get() /
                          float(runParams.operatorRuns
                          * runParams.metroStepsPerOperatorRun)).mean()

    if runParams.returnOperator:
        if not runParams.enableOperator:
            raise Exception("Can only return operator when "
                            "runParams.enableOperator = True")
        else:
            rawOpVector = pimcKernel.operatorValues.get()
            RKR.operatorMean = np.empty((len(runParams.operators),
                runParams.nbrOfWalkers))

            if runParams.enableParallelizePath:
                nbrOfThreadsPerWalker = runParams.N / (2 ** runParams.S)

                for j in range(len(runParams.operators)):
                    for i in range(runParams.nbrOfWalkers):
                        RKR.operatorMean[j, i] = (rawOpVector[j +
                            i * len(runParams.operators) * nbrOfThreadsPerWalker+
                           np.array(range(nbrOfThreadsPerWalker)) *
                           len(runParams.operators)].mean())
            else:
                for j in range(len(runParams.operators)):
                    for i in range(runParams.nbrOfWalkers):
                        RKR.operatorMean[j, i] = rawOpVector[j + i *
                                len(runParams.operators)].mean()

    if runParams.returnCorrelator:
        if not runParams.enableCorrelator:
            raise Exception("Can only return correlator when "
                            "runParams.enableCorrelator = True")
        correlatorValues = pimcKernel.correlatorValues.get()
        RKR.correlatorMean = (correlatorValues.mean(axis = 0) /
                       float(runParams.operatorRuns * runParams.N))
        RKR.correlatorStandardError = (correlatorValues.std(axis = 0) /
                        np.sqrt(pimcKernel.runParams.nbrOfWalkers))

    #Sum the bin counts from all individual threads.
    if runParams.returnBinCounts:
        if not runParams.enableBins:
            raise Exception("Can only return bins when "
                            "runParams.enableBins = True")
        totBinCount = (runParams.nbrOfWalkers * runParams.operatorRuns
                * runParams.N)
        binVol = ((float(runParams.xmax - runParams.xmin)
                 / float(runParams.binResolutionPerDOF)) **
                 pimcKernel.system.DOF)
        # Normalize probability density
        RKR.binCountNormed = (pimcKernel.binCounts.get().astype(np.float32)
                              / (binVol * totBinCount))

    #Return paths, useful when recompiling after burn in
    if runParams.returnPaths:
        RKR.paths = pimcKernel.paths.get()

    return RKR

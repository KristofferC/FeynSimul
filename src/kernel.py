# This file is part of FeynSimul.
#
# FeynSimul is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FeynSimul is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FeynSimul.  If not, see <http://www.gnu.org/licenses/>.

import sys
import numpy as np
from time import time
from collections import defaultdict
import os
import copy

import pyopencl as cl
import pyopencl.array
import pyopencl.clrandom
import pyopencl.clmath

def humanReadableSize(size):
    """
    Converts a size in bytes to human readable form.

    @type size: Number
    @param size: Size in bytes

    @rtype: string
    @return: The size in human readable form.
    """
    if size <= 0 or size >= 10e18:
        return "%.3g" % size + " B"
    else:
        p = int(np.log(size) / np.log(1024))
        names = ["", "ki", "Mi", "Gi", "Ti", "Pi", "Ei"]
        return "%.3g" % (size / 1024.0 ** p) + " " + names[p] + "B"


class PIMCKernel:
    def __init__(self,ka):
        """
        Loads the kernel to the GPU and readies it to be run.

        The kernel code from kernel.src is read and modified with the necessary
        defines and dictionary replacements. The binary program that will be
        sent to the GPU is built and memory allocations on the GPU are done.

        @type ka: L{kernelArgs} class
        @param self._system: An instance of kernelArgs describing what kind of 
        kernel to build.
        """
        self._enableBins = ka.enableBins
        self._enableOperator = ka.enableOperator
        self._enableCorrelator = ka.enableCorrelator
        self._enablePathShift = ka.enablePathShift
        self._enableBisection = ka.enableBisection
        self._enableSingleNodeMove = ka.enableSingleNodeMove
        self._enableParallelizePath = ka.enableParallelizePath
        self._enableGlobalPath = ka.enableGlobalPath
        self._enableGlobalOldPath = ka.enableGlobalOldPath
        self._nbrOfWalkers = ka.nbrOfWalkers
        self._nbrOfWalkersPerWorkGroup = ka.nbrOfWalkersPerWorkGroup
        self._N = ka.N 
        self._S = ka.S 
        self._beta = ka.beta
        self._operatorRuns = ka.operatorRuns
        self._metroStepsPerOperatorRun = ka.metroStepsPerOperatorRun
        self._system = copy.copy(ka.system)
        
        if self._enableOperator:
            self._operators = ka.operators

        if self._enableCorrelator:
            self._correlators = ka.correlators
        
        if self._enableBins:
            self._binResolutionPerDOF = ka.binResolutionPerDOF
            self._xMin = ka.xMin
            self._xMax = ka.xMax

        if not (self._nbrOfWalkers / self._nbrOfWalkersPerWorkGroup *
               self._nbrOfWalkersPerWorkGroup == self._nbrOfWalkers):
            raise Exception("The number of walkers must be a multiple "
                            "of the number of walkers per work group.")

        if self._nbrOfWalkers < self._nbrOfWalkersPerWorkGroup:
            raise Exception("The number of walkers per group must be less "
                            "than or equal to the total number of walkers")

        if not (self._N != 0 and ((self._N & (self._N - 1)) == 0)):
            raise Exception("Path node number, N, must be a power of 2.")

        if self._enableBisection:
            if 2 ** self._S > self._N:
                raise Exception("2^S must be less than or equal to N.")

        if self._enableParallelizePath:
            if self._enableSingleNodeMove:
                raise Exception("Single node move not possible when "
                                "enableParallelizePath is True.")
            if self._enablePathShift:
                raise Exception("Shift path not possible when "
                                "enableParallelizePath is True.")

            nbrOfThreadsPerWalker = self._N / (2 ** self._S)
            nbrOfWorkgroups = self._nbrOfWalkers / self._nbrOfWalkersPerWorkGroup

            self._localSize = (nbrOfThreadsPerWalker,
                                    self._nbrOfWalkersPerWorkGroup)
            self._globalSize = (nbrOfWorkgroups * nbrOfThreadsPerWalker,
                                     self._nbrOfWalkersPerWorkGroup)
            self._nbrOfThreads = self._nbrOfWalkers * nbrOfThreadsPerWalker
        else:
            self._localSize = None
            self._globalSize = (self._nbrOfWalkers,)
            self._nbrOfThreads = self._nbrOfWalkers

        defines = ""

        if self._enableBins:
            defines += "#define ENABLE_BINS\n"

        if self._enableOperator:
            defines += "#define ENABLE_OPERATOR\n"

        if self._enableCorrelator:
            defines += "#define ENABLE_CORRELATOR\n"

        if self._enablePathShift:
            defines += "#define ENABLE_PATH_SHIFT\n"

        if self._enableBisection:
            defines += "#define ENABLE_BISECTION\n"

        if self._enableSingleNodeMove:
            defines += "#define ENABLE_SINGLE_NODE_MOVE\n"

        if self._enableGlobalPath:
            defines += "#define ENABLE_GLOBAL_PATH\n"

        if self._enableGlobalOldPath:
            if not self._enableBisection:
                raise Exception("Old Path is not used when "
                                "bisection is not enabled")
            defines += "#define ENABLE_GLOBAL_OLD_PATH\n"

        if self._enableParallelizePath:
            defines += "#define ENABLE_PARALELLIZE_PATH\n"

        defines += "#define DOF_ARGUMENT_DECL "
        for i in range(self._system.DOF):
            defines += "float x" + str(i + 1)
            if (i != self._system.DOF - 1):
                defines += ","
        defines += "\n"

        defines += "#define DOF_ARGUMENT_DATA "
        for i in range(self._system.DOF):
            defines += "pathPointPtr[" + str(i) + " * " + str(self._N) + "]"
            if (i != self._system.DOF - 1):
                defines += ","
        defines += "\n"

        if self._enableOperator:
            operatorCode = ""
            nbrOfOperatorsZeros = ""
            for i in range(0, len(self._operators)):
                operatorCode += ("opAccum[" + str(i) + "] += "
                                 + self._operators[i] + ";\n")
                if i != 0:
                    nbrOfOperatorsZeros += ","
                nbrOfOperatorsZeros += "0.0f"

        if self._enableCorrelator:
            correlatorCode = ""
            nbrOfCorrelatorsOnes = ""
            for i in range(0, len(self._correlators)):
                correlatorCode += ("cokarod[" + str(i) + "] *= "
                                   + self._correlators[i] + ";\n")
                if i != 0:
                    nbrOfCorrelatorsOnes += ","
                nbrOfCorrelatorsOnes += "1.0f"

        class DictWithDefault(defaultdict):
            def __missing__(self, key):
                return key + str(" is not defined")

        # Replacements that will be done in the kernel code file will be stored
        # in this dictionary.
        replacements = DictWithDefault()

        replacements['nbrOfWalkersPerWorkGroup'] = self._nbrOfWalkersPerWorkGroup
        replacements['N'] = '%d' % self._N
        replacements['potential'] = self._system.potential
        replacements['epsilon'] = '%ef' % (self._beta / float(self._N))
        replacements['epsilon_inv2'] = '%ef' % ((float(self._N) / self._beta) ** 2)
        replacements['pathSize'] = '%d' % (self._system.DOF * self._N)
        replacements['defines'] = defines
        replacements['userCode'] = self._system.userCode
        replacements['DOF'] = '%d' % self._system.DOF
        replacements['operatorRuns'] = '%d' % self._operatorRuns
        replacements['metroStepsPerOperatorRun'] = ('%d'
                % self._metroStepsPerOperatorRun)

        if self._enableOperator:
            replacements['operatorCode'] = operatorCode
            replacements['nbrOfOperators'] = '%d' % len(self._operators)
            replacements['nbrOfOperatorsZeros'] = nbrOfOperatorsZeros
            replacements['opNorm'] = '%ef' % (1.0 / float(self._operatorRuns 
            * self._N * self._nbrOfWalkers / self._nbrOfThreads))

        if self._enableCorrelator:
            replacements['correlatorCode'] = correlatorCode
            replacements['nbrOfCorrelators'] = '%d' % len(self._correlators)
            replacements['nbrOfCorrelatorsOnes'] = nbrOfCorrelatorsOnes

        if self._enableBisection:
            replacements['PI'] = '%ef' % 3.141592653589793
            replacements['2_POW_S'] = '%d' % 2 ** self._S
            replacements['S'] = '%d' % self._S

        if self._enableSingleNodeMove:
            replacements['alpha'] = '%ef' % self._alpha

        if self._enablePathShift:
            replacements['PSAlpha'] = '%ef' % self._PSAlpha

        if self._enableBins:
            replacements['xMin'] = '%ef' % self._xMin
            replacements['xMax'] = '%ef' % self._xMax
            replacements['binsPerPart'] = '%d' % self._binResolutionPerDOF
            replacements['nbrOfBins'] = '%d' % self._binResolutionPerDOF ** self._system.DOF
            replacements['invBinSize'] = '%ef' % (float(self._binResolutionPerDOF) /
                                         float(self._xMax - self._xMin))

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
        self._ctx = cl.create_some_context()
        queueProperties = cl.command_queue_properties.PROFILING_ENABLE
        self._queue = cl.CommandQueue(self._ctx,
                                           properties=queueProperties)

        #Build the program and identify metropolis as the kernel
        self._prg = (cl.Program(self._ctx, kernelCode)
                         .build(options="-cl-fast-relaxed-math"))
        self._kernel = self._prg.metropolis

        #Initial paths are created (the initial path vector is filled with zeros,
        #meaning no movement of the particles)
        try:
            self._paths = cl.array.zeros(self._queue,
                              (self._nbrOfWalkers, self._N * self._system.DOF),
                              np.float32)

            #Buffer for storing number of accepted values and
            #seeds for the xorshfitPRNG
            self._accepts = cl.array.zeros(self._queue,
                    (self._nbrOfThreads, ), np.uint32)

            #np.random.seed(0)
            self._seeds = cl.array.to_device(self._queue,
                             (np.random.randint(0, high = 2 ** 31 - 1,
                              size = (self._nbrOfThreads + 1, 4))
                              ).astype(np.uint32))
            if self._enableOperator:
                #pyopencl.array objects are created for storing
                #the calculated operator means from each thread
                self._operatorValues = cl.array.zeros(self._queue,
                        self._nbrOfThreads * len(self._operators), np.float32)

            if self._enableCorrelator:
                #pyopencl.array objects are created for storing
                #the calculated operator means from each thread
                self._correlatorValues = cl.array.zeros(self._queue,
                        (self._nbrOfWalkers, len(self._correlators), self._N / 2),
                        np.float32)

            if self._enableGlobalOldPath:
                self._oldPath = cl.array.zeros(self._queue,
                                    (self._nbrOfThreads,
                                    (2 ** self._S - 1) * self._system.DOF), np.float32)

            #A buffer for the multidimensional array for storing
            #the resulting number of bin counts
            if self._enableBins:
                binTouple = (self._binResolutionPerDOF,)
                for i in range(1, self._system.DOF):
                    binTouple += (self._binResolutionPerDOF,)
                self._binCounts = cl.array.zeros(self._queue,
                        binTouple, np.uint32)

        except pyopencl.MemoryError:
            raise Exception("Unable to allocate global "
                            "memory on device, out of memory?")

    def run(self):
        """
        Runs the kernel with the self._system and the run parameters.

        This will start the kernel on the GPU. When the GPU is done desired
        data is fetched.

        @type pimcKernel: L{pimc_kernel.PIMCKernel}
        @param pimcKernel: The class that is returned from L{loadKernel}.
        @rtype: L{RunKernelResults}
        @return: Contains the results of running the kernel.
        """
        #Run kernel. The buffers for all pyopencl.array objects are passed
        #by the .data parameter.

        args = [self._paths.data,
                self._accepts.data,
                self._seeds.data]

        if self._enableGlobalOldPath:
            args += [self._oldPath.data]

        if self._enableOperator:
            args += [self._operatorValues.data]

        if self._enableCorrelator:
            args += [self._correlatorValues.data]

        if self._enableBins:
            args += [self._binCounts.data]

        #print(self._globalSize)#(nbrOfWorkgroups* nbrOfThreadsPerWalker, RP.nbrOfWalkersPerWorkGroup)
        #print(self._localSize)#(nbrOfThreadsPerWalker,RP.nbrOfWalkersPerWorkGroup)
        #print(self._getStats())

        startTime = time()
        self._kernelObj = self._kernel(self._queue,
                                          self._globalSize,
                                          self._localSize,
                                          *args)

        #Wait until the threads have finished and then calculate total run time
        try:
            self._kernelObj.wait()
        except pyopencl.RuntimeError:
            if time() - startTime > 5.0:
                raise Exception("Kernel runtime error. Over 5 seconds had passed "
                                "when kernel aborted.")
            else:
                raise

    def getPaths(self):
        return self._paths.get()
        
    def getOperators(self):
        if not self._enableOperator:
                raise Exception("Can only return operator when "
                                "enableOperator = True")
        rawOpVector = self._operatorValues.get()
        operatorMean = np.empty((len(self._operators),
            self._nbrOfWalkers))

        if self._enableParallelizePath:
            nbrOfThreadsPerWalker = self._N / (2 ** self._S)

            for j in range(len(self._operators)):
                for i in range(self._nbrOfWalkers):
                    operatorMean[j, i] = (rawOpVector[j +
                    i * len(self._operators) * nbrOfThreadsPerWalker +
                    np.array(range(nbrOfThreadsPerWalker)) *
                    len(self._operators)].mean())
        else:
            for j in range(len(self._operators)):
                for i in range(self._nbrOfWalkers):
                    operatorMean[j, i] = rawOpVector[j + i *
                        len(self._operators)].mean()
        return operatorMean

    def getCorrelator(self):
        if not self._enableCorrelator:
            raise Exception("Can only return correlator when "
                            "kernalArg.enableCorrelator = True")
        correlatorValues = self._correlatorValues.get()
        correlatorMean = (correlatorValues.mean(axis = 0) /
                       float(operatorRuns * N))
        correlatorStandardError = (correlatorValues.std(axis = 0) /
                        np.sqrt(self._nbrOfWalkers))
        return (correlatorMean, correlatorStandardError)

    def getBinCounts(self):
        if not self._enableBins:
            raise Exception("Can only return bins when "
                            "kernalArg.enableBins = True")
        #Sum the bin counts from all individual threads.
        totBinCount = (self._nbrOfWalkers * self._operatorRuns
                * N)
        binVol = ((float(self._xMax - self._xMin)
                 / float(self._binResolutionPerDOF)) **
                 self._system.DOF)
        # Normalize probability density
        return (self._binCounts.get().astype(np.float32)
                              / (binVol * totBinCount))
    def getRunTime(self):
        return 1e-9 * (self._kernelObj.profile.end - 
        self._kernelObj.profile.start)

    def getAcceptanceRate(self):
        return (self._accepts.get() / float(self._operatorRuns
        * self._metrostepsPerOperatorRun)).mean()

    def getStats(self):
        """
        Returns information about memory usage of the loaded PIMCKernel.
        Global memory usage is just estimated based on RunParameters and might
        thus not be completely accurate.

        @rtype:   string
        @return:  The information in human readable text format
        """

        ret = ""
        devices = self._prg.get_info(cl.program_info.DEVICES)
        if len(devices) != 1:
            raise Exception("Expected the number of " +
                            "devices to be 1, it was " + str(len(devices)))
        dev = devices[0]
        usedGlobalMemory = 0
        usedGlobalMemory += (self._nbrOfThreads + 1) * 4 * 4  # seeds
        usedGlobalMemory += (self._nbrOfThreads) * 4  # accepts
        usedGlobalMemory += (self._nbrOfWalkers * self._N *
                             self._self._system.DOF * 4)  # path
        usedGlobalMemory += (self._nbrOfThreads
                            * self._nbrOfOperators * 4)  # operator
        if self._enableGlobalOldPath:
            usedGlobalMemory += (self._nbrOfThreads *
                                (2 ** self._S - 1)
                                 * self._self._system.DOF * 4)
        if self._enableBins:
            usedGlobalMemory += (self._binsPerPart **
                                self._self._system.DOF * 4)

        ret += ("Global memory (used/max): " +
                humanReadableSize(usedGlobalMemory) + " / " +
                humanReadableSize(dev.get_info(cl.device_info.GLOBAL_MEM_SIZE))
                + "\n")

        ret += ("Local memory (used/max): " +
            humanReadableSize(self._prg.metropolis.get_work_group_info(
                                cl.kernel_work_group_info.LOCAL_MEM_SIZE, dev))
                                + " / " +
                                humanReadableSize(
                                dev.get_info(cl.device_info.LOCAL_MEM_SIZE))
                                + "\n")

        if self._enableParallelizePath:
            ret += ("Workgroup size (used/device max/kernel max): " +
                    str(self._N / (2 ** self._S) *
                    self._nbrOfWalkersPerWorkGroup)
                    + " / " +
                    str(dev.get_info(cl.device_info.MAX_WORK_GROUP_SIZE))
                    + " / " +
                    str(self._kernel.get_work_group_info(
                    cl.kernel_work_group_info.WORK_GROUP_SIZE, dev)) + "\n")
        ret += ("Workgroup dimensions (used/max): " +
                str(self._localSize) + " / " +
                str(dev.get_info(cl.device_info.MAX_WORK_ITEM_SIZES)) + "\n")
        ret += ("Number of workgroups (used): " +
                str(self._nbrOfWalkers
                    / self._nbrOfWalkersPerWorkGroup))
        return ret

    def exceedsLimits(self):
        """
        Returns true if information as from getStats indicates a violation
        of limits. Global memory usage is just estimated based on RunParameters
        and might thus not be completely accurate.

        @rtype:   bool
        @return:  True if limits are exceeded
        """
        devices = self._prg.get_info(cl.program_info.DEVICES)
        if len(devices) != 1:
            raise Exception("Expected the number of devices to be 1, it was "
                            + str(len(devices)))
        dev = devices[0]
        usedGlobalMemory = 0
        usedGlobalMemory += (self._nbrOfThreads + 1) * 4 * 4  # seeds
        usedGlobalMemory += self._nbrOfThreads * 4  # accepts
        usedGlobalMemory += (self._nbrOfWalkers *
                            self._N * self._self._system.DOF * 4)  # path
        usedGlobalMemory += (self._nbrOfThreads
                            * self._nbrOfOperators * 4)  # operator

        if self._enableGlobalOldPath:
            usedGlobalMemory += (self._nbrOfThreads *
                                 (2 ** self._S - 1)
                                 * self._self._system.DOF * 4)

        if self._enableBins:
            usedGlobalMemory += (self._binsPerPart
                                ** self._self._system.DOF * 4)

        if(usedGlobalMemory > dev.get_info(cl.device_info.GLOBAL_MEM_SIZE)):
            return True

        if(self._kernel.get_work_group_info(
           cl.kernel_work_group_info.LOCAL_MEM_SIZE, dev) >
           dev.get_info(cl.device_info.LOCAL_MEM_SIZE)):
            return True

        if self._enableParallelizePath:
            if ((self._N / (2 ** self._S) *
                self._nbrOfWalkersPerWorkGroup) >
                dev.get_info(cl.device_info.MAX_WORK_GROUP_SIZE)):
                return True
            if ((self._N / (2 ** self._S) *
                self._nbrOfWalkersPerWorkGroup) >
                self._kernel.get_work_group_info(
                cl.kernel_work_group_info.WORK_GROUP_SIZE, dev)):
                return True

        for i in range(len(self._localSize)):
            if (self._localSize[i] >
                dev.get_info(cl.device_info.MAX_WORK_ITEM_SIZES)[i]):
                return True

        return False

        #awaiting new version of pyopencl
'''def getBuildLog(self):

        devices=self._prg.get_info(cl.program_info.DEVICES)
        if len(devices)!=1:
            raise Exception("Expected the number of devices to be 1, it was "
+str(len(devices)))
        return self._prg.get_build_info(devices[0],
cl.program_build_info.LOG)'''

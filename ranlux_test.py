# Simple sandbox for RANLUX

import numpy as np
import pyopencl as cl
from time import time
from collections import defaultdict
import os
import copy

import pyopencl as cl
import pyopencl.array
import pyopencl.clrandom
import pyopencl.clmath

luxuaryFactor = 2
enableDouble = False
nbrOfWalkers = 448*32
N = 128
localSize = None
globalSize = (nbrOfWalkers,)
nbrOfThreads = nbrOfWalkers
randsPerThread = 400 # must be a multiple of 4!

defines = ""

programBuildOptions = "-cl-fast-relaxed-math"

if enableDouble:
    defines += "#define ENABLE_DOUBLE\n"
else:
    programBuildOptions += " -cl-single-precision-constant"

class DictWithDefault(defaultdict):
    def __missing__(self, key):
        return key + str(" is not defined")
        
replacements = DictWithDefault()
replacements['defines'] = defines
replacements['luxuaryFactor'] = str(luxuaryFactor)
replacements['randsPerThread'] = str(randsPerThread)

ctx = cl.create_some_context()
queueProperties = cl.command_queue_properties.PROFILING_ENABLE
queue = cl.CommandQueue(ctx, properties=queueProperties)
mf = cl.mem_flags

#initKernelCode_r = open(os.path.dirname(__file__) + 'ranlux_init_kernel.c', 'r').read()
#initKernelCode = initKernelCode_r % replacements
dummyBuffer = np.zeros(nbrOfThreads * 28, dtype=np.uint32)

ins = cl.array.to_device(queue, (np.random.randint(0, high = 2 ** 31 - 1, size = (nbrOfThreads))).astype(np.uint32))
ranluxcltab = cl.Buffer(ctx, mf.READ_WRITE, size=0, hostbuf=dummyBuffer)
#prg = (cl.Program(ctx, initKernelCode).build(options=programBuildOptions))
#kernel = prg.ranlux_init_kernel

#kernelObj = kernel(queue, globalSize, localSize, ins.data, )
#kernelObj.wait()

kernelCode_r = open(os.path.dirname(__file__) + 'ranlux_test_kernel.c', 'r').read()
kernelCode = kernelCode_r % replacements

prg = (cl.Program(ctx, kernelCode).build(options=programBuildOptions))

kernel_1 = prg.ranlux_init_kernel

kernelObj_1 = kernel_1(queue, globalSize, localSize, ins.data, ranluxcltab)
kernelObj_1.wait()

kernel_2 = prg.ranlux_test_kernel

randomsOut = cl.array.zeros(queue, nbrOfThreads * randsPerThread, np.float64 if enableDouble else np.float32)

kernelObj_2 = kernel_2(queue, globalSize, localSize, randomsOut.data, ranluxcltab)
kernelObj_2.wait()

resultingNumbers = randomsOut.get()


print '--- Running %d threads with %d random numbers per thread ---' % (nbrOfThreads, randsPerThread)
print 'Total number of rands: %d' % (nbrOfThreads * randsPerThread)
print 'Initiation time: %f' % (1e-9 * (kernelObj_1.profile.end - kernelObj_1.profile.start))
print 'Run time: %f' % (1e-9 * (kernelObj_2.profile.end - kernelObj_2.profile.start))
print 'RANLUX luxuary factor: %d' % luxuaryFactor
if enableDouble:
    print 'Float precision: double'
else:
    print 'Float precision: single'
print 'Mean: %f' % np.mean(resultingNumbers)
print 'Standard deviation: %f' % np.std(resultingNumbers)

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


# the histogram of the data
n, bins, patches = plt.hist(resultingNumbers, 50, normed=1)

plt.show()


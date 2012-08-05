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
enableRanlux = True
useRanluxInt = False
enableDouble = False
returnRandoms = True
enablePlot = False
printKernelCode = False
nbrOfWalkers = 448*32*4
localSize = None
globalSize = (nbrOfWalkers,)
nbrOfThreads = nbrOfWalkers
randsPerThread = 10000 # must be a multiple of 4!

defines = ""

programBuildOptions = "-cl-fast-relaxed-math"

if enableDouble:
    defines += "#define ENABLE_DOUBLE\n"
else:
    programBuildOptions += " -cl-single-precision-constant"

if returnRandoms:
    defines += "#define RETURN_RANDOMS\n"
    
if enableRanlux:
    defines += "#define ENABLE_RANLUX\n"
    
if useRanluxInt:
    defines += "#define USE_RANLUXCL_INT\n"

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

#initKernelCode_r = open(os.path.dirname(__file__) + 'ranlux_init_kernel.c', 'r').read()
#initKernelCode = initKernelCode_r % replacements

if enableRanlux:
    mf = cl.mem_flags
    dummyBuffer = np.zeros(nbrOfThreads * 28, dtype=np.uint32)
    ins = cl.array.to_device(queue, (np.random.randint(0, high = 2 ** 31 - 1, size = (nbrOfThreads))).astype(np.uint32))
    ranluxcltab = cl.Buffer(ctx, mf.READ_WRITE, size=0, hostbuf=dummyBuffer)
else:
    ins = cl.array.to_device(queue, (np.random.randint(0, high = 2 ** 31 - 1, size = (nbrOfThreads,4))).astype(np.uint32))
#prg = (cl.Program(ctx, initKernelCode).build(options=programBuildOptions))
#kernel = prg.ranlux_init_kernel

#kernelObj = kernel(queue, globalSize, localSize, ins.data, )
#kernelObj.wait()

kernelCode_r = open(os.path.dirname(__file__) + 'ranlux_test_kernel.c', 'r').read()
kernelCode = kernelCode_r % replacements

if printKernelCode:
    print '-'*20
    print kernelCode
    print '-'*20

prg = (cl.Program(ctx, kernelCode).build(options=programBuildOptions))

#kernel_1 = prg.ranlux_init_kernel

#kernelObj_1 = kernel_1(queue, globalSize, localSize, ins.data, ranluxcltab)
#kernelObj_1.wait()


if enableRanlux:
    kernel = prg.ranlux_test_kernel
else:
    kernel = prg.xorshift_test_kernel

if returnRandoms:
    if useRanluxInt and enableRanlux:
        randomsOut = cl.array.zeros(queue, nbrOfThreads * randsPerThread, np.uint32)
    else:
        randomsOut = cl.array.zeros(queue, nbrOfThreads * randsPerThread, np.float64 if enableDouble else np.float32)
        
        
    if enableRanlux:
        kernelObj = kernel(queue, globalSize, localSize, ins.data, randomsOut.data, ranluxcltab)
    else:
        kernelObj = kernel(queue, globalSize, localSize, ins.data, randomsOut.data)
else:
    if enableRanlux:
        kernelObj = kernel(queue, globalSize, localSize, ins.data, ranluxcltab)
    else:
        kernelObj = kernel(queue, globalSize, localSize, ins.data)

kernelObj.wait()

if returnRandoms:
    resultingNumbers = randomsOut.get()


print '--- Running %d threads with %d random numbers per thread ---' % (nbrOfThreads, randsPerThread)
print 'Total number of rands: %d' % (nbrOfThreads * randsPerThread)
if enableRanlux:
    print 'PRNG: RANLUX (luxuary factor %d)' % luxuaryFactor
else:
    print 'PRNG: xorshift'
#print 'Initiation time: %f' % (1e-9 * (kernelObj_1.profile.end - kernelObj_1.profile.start))
print 'Run time: %f' % (1e-9 * (kernelObj.profile.end - kernelObj.profile.start))
if useRanluxInt and enableRanlux:
    print 'Data type: uint 32bit'
else:
    if enableDouble:
        print 'Data type: float 64bit'
    else:
        print 'Data type: float 32bit'
if returnRandoms:
    print 'Mean: %f' % np.mean(resultingNumbers, dtype=np.float64)
    print 'Standard deviation: %f (expected from uniformly dist: 1/sqrt(12) = %f)' % (np.std(resultingNumbers, dtype=np.float64), (1.0/np.sqrt(12)))

if enablePlot and returnRandoms:
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt


    # the histogram of the data
    n, bins, patches = plt.hist(resultingNumbers, 50, normed=1)

    plt.show()


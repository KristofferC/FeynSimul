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
ranluxIntMax = 2**32-1
enableDouble = False
returnRandoms = True
enablePlot = False
printKernelCode = False
nbrOfWalkers = 8
ydim = 2
groupSize = 2
localSize = (groupSize,ydim)
globalSize = (nbrOfWalkers,ydim)
nbrOfThreads = nbrOfWalkers * ydim
randsPerThread = 1 # must be a multiple of 4!

defines = ""

programBuildOptions = "-cl-fast-relaxed-math -cl-nv-verbose"

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
replacements['ranluxIntMax'] = (str('%.8E' % ranluxIntMax) + 'f')

ctx = cl.create_some_context()
queueProperties = cl.command_queue_properties.PROFILING_ENABLE
queue = cl.CommandQueue(ctx, properties=queueProperties)

if enableRanlux:
    mf = cl.mem_flags
    dummyBuffer = cl.array.to_device(queue, np.zeros(nbrOfThreads * 28, dtype=np.uint32))
    ins = cl.array.to_device(queue, (np.random.randint(0, high = 2 ** 31 - 1, size = (nbrOfThreads/groupSize))).astype(np.uint32))
else:
    ins = cl.array.to_device(queue, (np.random.randint(0, high = 2 ** 31 - 1, size = (nbrOfThreads,4))).astype(np.uint32))

kernelCode_r = open(os.path.dirname(__file__) + 'ranlux_test_kernel.c', 'r').read()
kernelCode = kernelCode_r % replacements

if printKernelCode:
    print '-'*20
    print kernelCode
    print '-'*20

prg = (cl.Program(ctx, kernelCode).build(options=["-I", "."]))

kernel_1 = prg.ranlux_test_kernel_init
kernelObj_1 = kernel_1(queue, globalSize, localSize, ins.data, dummyBuffer.data)
kernelObj_1.wait()
queue.finish()
if enableRanlux:
    #kernel_init = prg.ranlux_test_kernel_init
    #kernelObj_init = kernel_init(queue, globalSize, localSize, ins.data, ranluxcltab)
    #kernelObj_init.wait()
    kernel = prg.ranlux_test_kernel
else:
    kernel = prg.xorshift_test_kernel
if returnRandoms:
    if useRanluxInt and enableRanlux:
        randomsOut = cl.array.zeros(queue, nbrOfThreads * randsPerThread, np.uint32)
    else:
        randomsOut = cl.array.zeros(queue, nbrOfThreads * randsPerThread, np.float64 if enableDouble else np.float32)
        randCountG = cl.array.zeros(queue, nbrOfThreads, np.uint32)
        
    if enableRanlux:
        testOut = cl.array.zeros(queue, nbrOfThreads * randsPerThread, np.uint32)
        #kernelObj = kernel(queue, globalSize, localSize, ins.data, randomsOut.data, testOut.data, dummyBuffer.data, randCountG.data)
    else:
        kernelObj = kernel(queue, globalSize, localSize, ins.data, randomsOut.data)
else:
    if enableRanlux:
        kernelObj = kernel(queue, globalSize, localSize, ins.data, dummyBuffer.data)
    else:
        kernelObj = kernel(queue, globalSize, localSize, ins.data)

#print dummyBuffer.get()
for i in range(0,8):
    """
    ins = cl.array.to_device(queue,ins.get(),np.uint32)
    randomsOut = cl.array.zeros(queue, nbrOfThreads * randsPerThread, dtype=np.uint32)
    testOut = cl.array.zeros(queue, nbrOfThreads * randsPerThread, dtype=np.uint32)
    dummyBuffer = cl.array.to_device(queue, dummyBuffer.get(), dtype=np.uint32)
    randCountG = cl.array.to_device(queue, randCountG.get(), dtype=np.uint32)
    """
    kernelObj = kernel(queue, globalSize, localSize, ins.data, randomsOut.data, testOut.data, dummyBuffer.data, randCountG.data)
    kernelObj.wait()
    queue.finish()
    if returnRandoms:
        resultingNumbers = randomsOut.get()
        resultingTest = testOut.get()
        if len(resultingTest) == len(resultingNumbers):
            for i in randCountG.get():
                print i, 
            print '\n',
            #for i in ins.get():
            #    print i, 
            #print '\n',
            for i in resultingTest:
                print '%d\t    ' % i, 
            print '\n',
            for i in resultingNumbers:
                if useRanluxInt:
                    print i,
                else:
                    print '%.5f' % i, 
        else:
            'Woopsie, something went wrong :('
    print '\n',
#print dummyBuffer.get()
#for i in range(0,len(resultingNumbers)):
#    if resultingNumbers[0] == resultingNumbers[i]:
#        print i

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
    if useRanluxInt and enableRanlux:
        print 'Mean: %f \t(normed: %f)' % (np.mean(resultingNumbers, dtype=np.float64), (np.mean(resultingNumbers, dtype=np.float64)/float(ranluxIntMax)))
        print 'Normed Standard deviation: %f (expected from uniformly dist: 1/sqrt(12) = %f)' % ((np.std(resultingNumbers, dtype=np.float64)/float(ranluxIntMax)), (1.0/np.sqrt(12)))
    else:
        print 'Mean: %.8f' % np.mean(resultingNumbers, dtype=np.float64)
        print 'Standard deviation: %f (expected from uniformly dist: 1/sqrt(12) = %f)' % (np.std(resultingNumbers, dtype=np.float64), (1.0/np.sqrt(12)))

if enablePlot and returnRandoms:
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt


    # the histogram of the data
    n, bins, patches = plt.hist(resultingNumbers, 50, normed=1)

    plt.show()


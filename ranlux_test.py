# Simple sandbox for RANLUX

import numpy as np
import pyopencl as cl

luxuaryFactor = 0
enableDouble = False
nbrOfWalkers = 448*32
N = 128
localSize = None
globalSize = nbrOfWalkers
nbrOfThreads = nbrOfWalkers
randsPerThread = 100

defines = ""

programBuildOptions = "-cl-fast-relaxed-math"

if enableDouble:
    defines += "#define ENABLE_DOUBLE\n"
else:
programBuildOptions += "-cl-single-precision-constant"

replacements['defines'] = defines
replacements['luxuaryFactor'] = luxuaryFactor
replacements['randsPerThread'] = randsPerThread

ctx = cl.create_some_context()
queueProperties = cl.command_queue_properties.PROFILING_ENABLE
queue = cl.CommandQueue(ctx, properties=queueProperties)

initKernelCode_r = open(os.path.dirname(__file__) + 'ranlux_init_kernel.c', 'r').read()
initKernelCode = initKernelCode_r % replacements

ins = cl.array.to_device(queue, (np.random.randint(0, high = 2 ** 31 - 1, size = (nbrOfThreads, 28))).astype(np.uint32))
prg = (cl.Program(ctx, initKernelCode).build(options=programBuildOptions))
kernel = prg.ranlux_test_kernel

kernelObj = kernel(queue, globalSize, localSize, ins.data)
kernelObj.wait()

kernelCode_r = open(os.path.dirname(__file__) + 'ranlux_test_kernel.c', 'r').read()
kernelCode = kernelCode_r % replacements

prg = (cl.Program(ctx, kernelCode).build(options=programBuildOptions))
kernel = prg.ranlux_test_kernel

ranluxcltab = cl.array.zeros(queue, nbrOfThreads * randsPerThread, np.float64 if enableDouble else np.float32)

kernelObj = kernel(queue, globalSize, localSize, ranluxcltab.data)
kernelObj.wait()


resultingNumbers = ranluxcltab.get()


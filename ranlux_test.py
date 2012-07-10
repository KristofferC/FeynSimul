# Simple sandbox for RANLUX

import numpy as np
import pyopencl as cl

def humanReadableSize(size):
    """
    Converts a size in bytes to human readable form.

    :param number size: Size in bytes to convert

    :rtype: string
    :returns: The size in human readable form.
    """
    if size <= 0 or size >= 10e18:
        return "%.3g" % size + " B"
    else:
        p = int(np.log(size) / np.log(1024))
        names = ["", "ki", "Mi", "Gi", "Ti", "Pi", "Ei"]
        return "%.3g" % (size / 1024.0 ** p) + " " + names[p] + "B"

luxuaryFactor = 0
enableDouble = False
nbrOfWalkers = 448
N = 128
localSize = None
globalSize = nbrOfWalkers
nbrOfThreads = nbrOfWalkers
numbersPerThread = 1000

defines = ""

if enableDouble:
    defines += "#define ENABLE_DOUBLE\n"

replacements['luxuaryFactor'] = luxuaryFactor
replacements['defines'] = defines

ctx = cl.create_some_context()
queueProperties = cl.command_queue_properties.PROFILING_ENABLE
queue = cl.CommandQueue(ctx, properties=queueProperties)

programBuildOptions = "-cl-fast-relaxed-math -cl-single-precision-constant"

kernelCode_r = open(os.path.dirname(__file__) + 'ranlux_test_kernel.c', 'r').read()
kernelCode = kernelCode_r % replacements

prg = (cl.Program(ctx, kernelCode).build(options=programBuildOptions))
kernel = prg.ranlux_test_kernel

randomNumbers = cl.array.zeros(queue, nbrOfThreads, np.float64 if enableDouble else np.float32)

kernelObj = kernel(queue, globalSize, localSize, randomNumbers.data)
kernelObj.wait()


resultingNumbers = randomNumbers.get()






ret = ""
devices = prg.get_info(cl.program_info.DEVICES)
if len(devices) != 1:
    raise Exception("Expected the number of " +
                    "devices to be 1, it was " + str(len(devices)))
dev = devices[0]
ret += ("Global memory (used/max): " +
        humanReadableSize(getGlobalMemory()) + " / " +
        humanReadableSize(dev.get_info(cl.device_info.GLOBAL_MEM_SIZE))
        + "\n")

ret += ("Local memory (used/max): " +
    humanReadableSize(prg.ranlux_test_kernel.get_work_group_info(
                        cl.kernel_work_group_info.LOCAL_MEM_SIZE, dev))
                        + " / " +
                        humanReadableSize(
                        dev.get_info(cl.device_info.LOCAL_MEM_SIZE))
                        + "\n")
print ret


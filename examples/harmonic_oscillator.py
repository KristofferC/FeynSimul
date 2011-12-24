# A simple run on a Harmonic Oscillator. The kernel is called once and the
# mean of the energy operator is printed. The bisection method with threads
# working in parallel is used

import os
import sys
import numpy as np

sys.path.append(sys.path[0] + "/../src/")
from run_params import *
from host import *

sys.path.append(sys.path[0] + "/../src/physical_systems/")
# Import the harmonic oscillator class
from harm_osc import *

# Create the 
HO_syst = HarmOsc()

HO_runparams = RunParams()
HO_runparams.nbrOfWalkers = 448*16
HO_runparams.N = 128
HO_runparams.S = 5
HO_runparams.beta = 15
HO_runparams.enableBisection = True
HO_runparams.operatorRuns = 100
HO_runparams.metroStepsPerOperatorRun = 10
HO_runparams.enableGlobalOldPath = False
HO_runparams.enableParallelizePath = True

operators = (HO_syst.energyOp,)

kernelEnv = loadKernel(HO_syst, HO_runparams, operators)
runKernelResults = runKernel(kernelEnv)

print runKernelResults.operatorMean


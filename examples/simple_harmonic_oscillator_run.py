# A simple run on a Harmonic Oscillator. The kernel is called once and the
# mean of the energy operator is printed. The bisection method with threads
# working in parallel is used. No thermalization is being done.

import os
import sys
import numpy as np

import pylab as pl

sys.path.append(sys.path[0] + "/../src/")
from run_params import *
from host import *

# Import the harmonic oscillator class
sys.path.append(sys.path[0] + "/../src/physical_systems/")
from harm_osc import *

# Create the physical system
HO_syst = HarmOsc()

# Set the run parameters
HO_runparams = RunParams()
HO_runparams.nbrOfWalkers = 448*2
HO_runparams.N = 256
HO_runparams.S = 6
HO_runparams.beta = 11
HO_runparams.operatorRuns = 20
HO_runparams.metroStepsPerOperatorRun = 20
HO_runparams.enableBisection = True
HO_runparams.enableParallelizePath = True
HO_runparams.binsEnabled = True
HO_runparams.returnBinCounts = True
HO_runparams.xmin = -4.0
HO_runparams.xmax = 4.0
HO_runparams.binResolutionPerDOF = 60


# Set the operator
operators = (HO_syst.energyOp,)

# Load kernel
HO_kernelEnvironment = loadKernel(HO_syst, HO_runparams, operators)

# Run kernel
HO_kernelResults = runKernel(HO_kernelEnvironment)

# Print the results
print "Mean: " + str(HO_kernelResults.operatorMean[0])
print "Standard error: " + str(HO_kernelResults.operatorStandardError)
print "Acceptance rate: " + str(HO_kernelResults.acceptanceRate)
print "Time for GPU to finish calculations: " + str(HO_kernelResults.runTime)

#Plot resulting wavefunction

#if plotWaveFunction and RP.binsEnabled:
  
# x interval for plot, + 0.5 * binSize to have each value in the middle of bins
binSize = (HO_runparams.xmax - HO_runparams.xmin) / HO_runparams.binResolutionPerDOF
x = np.linspace(HO_runparams.xmin, HO_runparams.xmax - binSize,
                HO_runparams.binResolutionPerDOF) + 0.5 * binSize
pl.plot(x, HO_kernelResults.binCountNormed,'*')
pl.show()


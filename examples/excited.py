# A simple run on a Harmonic Oscillator. The kernel is called once and the
# mean of the energy operator is printed. The bisection method with threads
# working in parallel is used. No thermalization is being done.

import os
import sys
import numpy as np
from sys import stdout
from time import sleep

import pylab as pl
from time import time
import numpy as np

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
HO_runparams.operatorRuns = 20
HO_runparams.metroStepsPerOperatorRun = 20
HO_runparams.enableBisection = True
HO_runparams.enableParallelizePath = True


# Set the operator
operators = (HO_syst.energyOp,)

energies = []
for beta in np.linspace(0.5,5,30):
    HO_runparams.beta = beta
    # Load kernel

    HO_kernelEnvironment = loadKernel(HO_syst, HO_runparams, operators)

###############################
## Burn in    
###############################
    startTime = time()
    burninTime = 3
    currTime = 0

    while currTime - startTime < burninTime:
        currTime = time()

        # Run kernel
        runKernel(HO_kernelEnvironment)
        k = currTime - startTime
        stdout.write("\r                                                             ")
        stdout.write ("\r Burning: " + "%2.3f" % k + " of " + str(burninTime))
        stdout.flush()
##################################

###############################
## Do run  
###############################
    stdout.write("\n")
    energyStore = 0.0
    runs = 0.0
    startTime = time()
    runTime = 5
    currTime = 0
    while currTime - startTime < runTime :
        currTime = time()
        runs += 1.0

        # Run kernel
        RKR = runKernel(HO_kernelEnvironment)
        energyStore += RKR.operatorMean[0]
        k = currTime - startTime
        stdout.write("\r                                                             ")
        stdout.write ("\r Running: " + "%2.3f" % k + " of " + str(runTime))
        stdout.flush()
        
    print "\n"
    print energyStore/runs
    energies.append(energyStore/runs)

print energies


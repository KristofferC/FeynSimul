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


# A simple run on a Harmonic Oscillator. The kernel is called once and the
# mean of the energy operator is printed. The bisection method with threads
# working in parallel is used. No thermalization is being done.

import os
import sys
import numpy as np

import pylab as pl

sys.path.append(sys.path[0] + "/../src/")
from kernel_args import *
from kernel import *

# Import the harmonic oscillator class
sys.path.append(sys.path[0] + "/../src/physical_systems/")
from harm_osc import *

# Set the run parameters
ka = KernelArgs()
ka.system = HarmOsc()
ka.nbrOfWalkers = 448*2
ka.N = 256
ka.S = 6
ka.beta = 11
ka.operatorRuns = 300
ka.enableOperator = True
ka.enableCorrelator = False
ka.metroStepsPerOperatorRun = 40
ka.enableBisection = True
ka.enablePathShift = False
ka.enableSingleNodeMove = False
ka.enableParallelizePath = True
ka.enableGlobalPath = False
ka.enableGlobalOldPath = False
ka.enableBins = True
ka.xMin = -4.0
ka.xMax = 4.0
ka.binResolutionPerDOF = 60
ka.nbrOfWalkersPerWorkGroup = 4
plotWaveFunction = True

# Set the operator
ka.operators  = (ka.system.energyOp,)

# Load kernel
kernel = PIMCKernel(ka)

# Run kernel
kernel.run()

# Print the results
print "Mean: " + str(np.mean(kernel.getOperators()))
#print "Standard error: " + str(HO_kernelResults.operatorStandardError)
#print "Acceptance rate: " + str(HO_kernelResults.acceptanceRate)
#print "Time for GPU to finish calculations: " + str(HO_kernelResults.runTime)

#Plot resulting wavefunction

if plotWaveFunction and ka.enableBins:
  
# x interval for plot, + 0.5 * binSize to have each value in the middle of bins
    binSize = (ka.xMax - ka.xMin) / ka.binResolutionPerDOF
    x = np.linspace(ka.xMin, ka.xMax - binSize,
                ka.binResolutionPerDOF) + 0.5 * binSize
    pl.plot(x,kernel.getBinCounts(),'*')
    pl.plot(x, 1/np.sqrt(np.pi) * np.exp(-x ** 2))
    pl.show()


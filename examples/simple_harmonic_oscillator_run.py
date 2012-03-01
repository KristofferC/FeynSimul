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


import numpy as np
import pylab as pl

from FeynSimul.kernel_args import *
from FeynSimul.kernel import *

# Import the harmonic oscillator class
from FeynSimul.physical_systems.harm_osc import *

# Set the run parameters
system = HarmOsc()
ka = KernelArgs(system = system,
                nbrOfWalkers = 448,
                N = 128,
                beta = 11.0,
                S = 6,
                enableOperator = True,
                enableCorrelator = False,
                metroStepsPerOperatorRun = 40,
                enableBisection = True,
                enablePathShift = False,
                enableSingleNodeMove = False,
                enableParallelizePath = True,
                enableGlobalPath = False,
                enableGlobalOldPath = False,
                enableBins = True,
                xMin = -3.5,
                xMax = 3.5,
                operatorRuns = 300,
                binResolutionPerDOF = 30,
                nbrOfWalkersPerWorkGroup = 1,
                operators = (system.energyOp,),
                correlators = ("x1",))

plotWaveFunction = True

# Load kernel
kernel = PIMCKernel(ka)

# Run kernel
kernel.run()

# Print the results
print "Ground state at: " + str(np.mean(kernel.getOperators()))


#Plot resulting wavefunction
if plotWaveFunction and ka.enableBins:

    binSize = (ka.xMax - ka.xMin) / ka.binResolutionPerDOF

    # x interval for plot, + 0.5 * binSize to have each value in the middle of bins
    x = np.linspace(ka.xMin, ka.xMax - binSize,
                ka.binResolutionPerDOF) + 0.5 * binSize
    pl.plot(x,kernel.getBinCounts(),'*', label="Simulated")

    # Analytical
    pl.plot(x, 1/np.sqrt(np.pi) * np.exp(-x ** 2), label="Analytical")
    pl.legend(loc='best')
    #pl.xlabel("x")
    #pl.title("$| \psi_0(x)|^2$ for Simple Harmonic Oscillator")
    pl.show()



#corrs = kernel.getCorrelator()[0]
#logDerCorr = -np.gradient(np.log(corrs[0]), ka.beta / ka.N)
#pl.plot(logDerCorr, '*')
#pl.show()

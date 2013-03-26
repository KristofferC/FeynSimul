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


# Does a run with the lm2m2 potential for 3 particles.

import os
import sys
import math
import numpy as np
from datetime import datetime
from time import time
from time import sleep
import csv

from FeynSimul.kernel_args import *
from FeynSimul.kernel import *

from FeynSimul.pimc_utils import *

# Import the harmonic oscillator class
from FeynSimul.physical_systems.lm2m2_cluster import *
enableDoubleOuter=False
system = Lm2m2_cluster(3,10e-3,verbose=False,enableDouble=enableDoubleOuter)
ka = KernelArgs(system = system,
                nbrOfWalkers = 64,
                N = 8,
                beta = system.beta,
                S = 6,
                enableOperator = True,
                enableCorrelator = False,
                metroStepsPerOperatorRun = 40,
                enableBisection = True,
                enablePathShift = False,
                enableSingleNodeMove = False,
                enableParallelizePath = True,
                enableGlobalPath = True,
                enableGlobalOldPath = True,
                enableBins = False,
                enableDouble = enableDoubleOuter,
                enableRanlux = False,
                luxuaryFactor = 2,
                ranluxIntMax = 2**32-1,
                operatorRuns = 100,
                nbrOfWalkersPerWorkGroup = 4,
                operators = (system.energyOp, system.meanSquaredRadiusOp))
# Time to run simul
#endTime = 60 * 60 * 24 * 14
endTime = 0 # Stop when final N has reached 'runsThisN' amount of runs

# How often to save paths.
savePathsInterval = 3000
ka.operators = (ka.system.energyOp, ka.system.meanSquaredRadiusOp)


def opRunsFormula(N, S):
    return max(2 ** 10 / 2 ** S, 1)

def mStepsPerOPRun(N, S):
    return 10

def runsPerN(N, S):
    return max(N / 8, 10)

startXList=np.random.uniform(size=(ka.nbrOfWalkers,ka.system.DOF)
        ,low=-1.0,high=1.0)*0.1*1e-10/ka.system.rUnit*np.sqrt(ka.system.m)*0

# Run the simulation function
modN(ka, startXList, savePathsInterval, "lm2m2_cluster", opRunsFormula
     , mStepsPerOPRun,  runsPerN, 512, simTime=endTime, finalN=1024*64,
     verbosity=2, continueRun=False, cont_S=12, cont_N=1024*16, cont_path='results/lm2m2_cluster/2013-02-20/13.51.01/pathsN16384episode1')

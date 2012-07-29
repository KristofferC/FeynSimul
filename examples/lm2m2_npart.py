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

sys.path.append(sys.path[0] + "/../src/")
from pimc_utils import *

sys.path.append(sys.path[0] + "/../src/physical_systems/")
import lm2m2_cluster

ka= KernelArgs()
ka.nbrOfWalkers = 64
ka.N = 8
ka.enablePathShift = False
ka.enableBisection = True
ka.enableSingleNodeMove = False
ka.enableGlobalPath = True
ka.enableGlobalOldPath = True
ka.enableParallelizePath = True
ka.enableBins = False
ka.nbrOfWalkersPerWorkGroup = 4
ka.operatorRuns = 100
ka.enableOperator = True
ka.enableCorrelator = False
ka.metroStepsPerOperatorRun = 40
ka.system = lm2m2_cluster.Lm2m2_cluster(3,15e-3,verbose=False)
ka.beta=ka.system.beta
# Time to run simul
endTime = 60 * 60 * 24 * 14

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
     , mStepsPerOPRun,  runsPerN, 512, runTime=endTime, finalN=1024*64,
     verbosity=2)

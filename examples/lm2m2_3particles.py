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

from FeynSimul.pimc_utils import *

from FeynSimul.physical_systems.lm2m2_3part import *

ka = KernelArgs(nbrOfWalkers = 64,
                N = 8,
                enablePathShift = False,
                enableBisection = True,
                enableSingleNodeMove = False,
                enableGlobalPath = True,
                enableGlobalOldPath = True,
                enableParallelizePath = True,
                enableBins = False,
                beta = 1000.0,
                nbrOfWalkersPerWorkGroup = 4,
                operatorRuns = 100,
                enableOperator = True,
                enableCorrelator = False,
                metroStepsPerOperatorRun = 40,
                system = Lm2m2_3part())

# Time to run simul
endTime = 60 * 60 * 24 * 14

# How often to save paths.
savePathsInterval = 3000
ka.operators = (ka.system.energyOp, ka.system.meanSquaredRadiusOp
              , ka.system.meanRadiusOp)


def opRunsFormula(N, S):
    return max(2 ** 10 / 2 ** S, 1)

def mStepsPerOPRun(N, S):
    return 10

def runsPerN(N, S):
    return max(N / 8, 10)

startXList=np.random.uniform(size=(ka.nbrOfWalkers,ka.system.DOF)
        ,low=-1.0,high=1.0)*0.01

# Run the simulation function
modN(ka, startXList, savePathsInterval, "lm2m2_3part", opRunsFormula
     , mStepsPerOPRun,  runsPerN, 512, runTime=endTime, finalN=1024*64,
     verbose=True)

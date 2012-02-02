# This file is part of FeynSimul.
#
# FeynSimul is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FeynSimul is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PUka.SE.  See the
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
from run_params import *
from host import *
from pimc_utils import *

sys.path.append(sys.path[0] + "/../src/physical_systems/")
from lm2m2_3part import *

ka.= RunParams()
ka.nbrOfWalkers = 64
ka.N = 1024 * 64
ka.getOperator = True
ka.enablePathShift = False
ka.enableBisection = True
ka.enableSingleNodeMove = False
ka.enableGlobalPath = True
ka.enableGlobalOldPath = True
ka.enableParallelizePath = True
ka.returnBinCounts = False
ka.beta = 1000
ka.nbrOfWalkersPerWorkGroup = 4


# Time to run simul
endTime = 60 * 60 * 24 * 14

# How often to save paths.
savePathsInterval = 3000
systemClass = Lm2m2_3part()
ka.operators = (systemClass.energyOp, systemClass.meanSquaredRadiusOp
              , systemClass.meanRadiusOp)


def opRunsFormula(N, S):
    return max(2 ** 10 / 2 ** S, 1)

def mStepsPerOPRun(N, S):
    return 10

def runsPerN(N, S):
    return max(N / 8, 10)

startXList=np.random.uniform(size=(ka.nbrOfWalkers,systemClass.DOF),low=-1.0,high=1.0)*0.01

# Run the simulation function
modN(ka, startXList, savePathsInterval, systemClass, endTime, "lm2m2_3part", opRunsFormula
     , mStepsPerOPRun, 8, runsPerN, 512)

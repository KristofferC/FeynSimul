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
from run_params import *
from host import *
from pimc_utils import *

sys.path.append(sys.path[0] + "/../src/physical_systems/")
from lm2m2_3part import *

RP = RunParams()
RP.nbrOfWalkers = 64
RP.N = 1024 * 64
RP.getOperator = True
RP.enablePathShift = False
RP.enableBisection = True
RP.enableSingleNodeMove = False
RP.enableGlobalPath = True
RP.enableGlobalOldPath = True
RP.enableParallelizePath = True
RP.returnBinCounts = False
RP.beta = 1000
RP.nbrOfWalkersPerWorkGroup = 4


# Time to run simul
endTime = 60 * 60 * 24 * 14

# How often to save paths.
savePathsInterval = 3000
systemClass = Lm2m2_3part()
RP.operators = (systemClass.energyOp, systemClass.meanSquaredRadiusOp
              , systemClass.meanRadiusOp)


def opRunsFormula(N, S):
    return max(2 ** 10 / 2 ** S, 1)

def mStepsPerOPRun(N, S):
    return 10

def runsPerN(N, S):
    return max(RP.N / 8, 10)

startXList=np.random.uniform(size=(RP.nbrOfWalkers,systemClass.DOF),low=-1.0,high=1.0)*0.01

# Run the simulation function
modN(RP, startXList, savePathsInterval, systemClass, endTime, "lm2m2_3part", opRunsFormula
     , mStepsPerOPRun, 8, runsPerN, 512)

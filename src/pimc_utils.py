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
# along with FeynSimul. If not, see <http://www.gnu.org/licenses/>.

import os
from datetime import datetime
from time import time
import sys
import copy

import csv
import numpy as np

from kernel_args import *
from kernel import *

def modN(ka, startXList, savePathsInterval, experimentName, opRunsFormula, 
        mStepsPerOPRun, runsPerN, maxWGSize, continueRun=False, finalN=-1,
        runTime=-1):

    startClock = time()
    pathChanges = 0
    timestamp  = datetime.now().strftime("%Y-%m-%d/%H.%M.%S")
    filename = "results/" + experimentName + "/" + timestamp + "/data"

    # Create directory to store results and paths
    if not os.path.exists("results/" + experimentName + "/" + timestamp):
        os.makedirs("results/" + experimentName + "/" + timestamp)
    
    kaMod = copy.copy(ka)

    if continueRun:
        kaMod.N = continueRun[1]
        kaMod.S = continueRun[2]
    else:
        kaMod.N = ka.N
        kaMod.S = 1
    
    firstNRun = True

    while finalN==-1 or kaMod.N <= finalN:
        ka.operatorRuns = opRunsFormula(kaMod.N, kaMod.S)
        ka.metroStepsPerOperatorRun = mStepsPerOPRun(kaMod.N, kaMod.S)

        kernel=PIMCKernel(kaMod)

        # If not first run create a new path by interpolating
        if firstNRun:
            if not continueRun:
                initialPaths = np.zeros((ka.nbrOfWalkers, kaMod.N * kaMod.system.DOF))
                for i in range(ka.nbrOfWalkers):
                    for j in range(kaMod.system.DOF):
                        initialPaths[i][j*kaMod.N:(j+1)*kaMod.N] = (np.ones(kaMod.N) *
                            startXList[i][j])
            else:
                initialPaths = np.zeros((ka.nbrOfWalkers, kaMod.N * kaMod.system.DOF))
                reader = csv.reader(open(continueRun[0]), delimiter='\t')
                k = 0
                for row in reader:
                    initialPaths[k] = map(float,row)
                    k += 1
            kernel.setPaths(initialPaths)
        else:
            newPaths = np.zeros((ka.nbrOfWalkers, kaMod.N * kaMod.system.DOF))
            for walker in range(ka.nbrOfWalkers):
                for DOF in range(kaMod.system.DOF):
                    secondIndices = DOF * kaMod.N + np.array(range(0, kaMod.N, 2))
                    # Nodes from the old path is copied to every second node in
                    # new path.
                    newPaths[walker, secondIndices] = oldPaths[walker
                            , secondIndices / 2]
                    # Linear interpolation to create new nodes in between nodes
                    # from old path.
                    # Fix periodic boundary conditions
                    rightIndex = secondIndices + 2
                    rightIndex[-1] -= kaMod.N
                    newPaths[walker, secondIndices + 1] = (newPaths[walker,
                        secondIndices] + newPaths[walker, rightIndex]) / 2.0
            # KE.paths.data.release()
            kernel.setPaths(newPaths)

        nRuns = 1
        runsThisN = runsPerN(kaMod.N, kaMod.S)
        while nRuns <= runsThisN or kaMod.N == finalN:
            if runTime != -1 and time() - startClock > runTime:
                return

            # Save paths
            if nRuns % savePathsInterval == 0 or nRuns == runsThisN :
                with open("results/" + experimentName + "/" + timestamp +
                         "/pathsN" + str(kaMod.N) + "episode" +
                         str(int(nRuns / savePathsInterval)), 'wb') as f:
                    csvWriter = csv.writer(f, delimiter='\t')
                    for aPath in kernel.getPaths():
                        csvWriter.writerow(aPath)
                print("Paths saved!")

            kernel.run()
            nRuns += 1
            pathChanges += kernel.getMetroStepsPerRun()
            output(filename, kaMod.N, time()-startClock, pathChanges
                   , kernel.getAcceptanceRate(), kernel.getOperators
                   , ka.beta, kaMod.S)

        #do preparations for next run that are not to be done first run
        if kaMod.N != finalN:
            oldPaths = kernel.getPaths()
            firstNRun = False
            # Change S to move acceptance rate in the right direction, some magic
            # numbers here.
            ar=kernel.getAcceptanceRate()
            if ar > 0.5:
                kaMod.S = min(i - 1, kaMod.S + 2)
            if 0.2 < ar < 0.5:
                kaMod.S = min(i - 1, kaMod.S + 1)
            if ar < 0.1:
                kaMod.S = max(1, kaMod.S - 1)
             # Make sure S is large enough to not use too many walkers in a WG
            while ka.nbrOfWalkersPerWorkGroup * kaMod.N / (2 ** kaMod.S) > maxWGSize:
                kaMod.S += 1
            kaMod.N *= 2

def output(filename, N, t, pathChanges, acceptanceRate
           , operatorMean, beta, S):

    print("N: " + str(N) + "\tS: " + str(S) + "\tbeta: " +
          str(beta)+"\tAR: " + str(acceptanceRate) + "\tOP:s " +
          str(np.mean(operatorMean, axis = 1)))

    with open(filename + ".tsv", 'a') as f_data:
        # Add a heading describing the columns if file is empty.
        if  os.path.getsize(filename + ".tsv") == 0:
            f_data.write("#Beta: " + str(beta) + "\n")
            f_data.write("#N\tTime\tpathChanges\tAR\tS")
            for i in range(len(operatorMean)):
                for j in range(len(operatorMean[0])):
                 f_data.write("\tOperator " + str(i) + ", Thread " + str(j) )
            f_data.write("\n")

        f_data.write(str(N)+"\t")
        f_data.write(str(t)+"\t")
        f_data.write(str(pathChanges)+"\t")
        f_data.write(str(acceptanceRate)+"\t")
        f_data.write(str(S)+"\t")
        for walkerOperators in operatorMean:
            for op in walkerOperators:
                f_data.write(str(op)+"\t")
        f_data.write("\n")

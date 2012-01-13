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


import os
from datetime import datetime
from time import time
import sys

import matplotlib.pyplot as plt
import csv
import numpy as np

sys.path.append(sys.path[0] + "/physical_systems/")
from lm2m2_4part import *
sys = Lm2m2_4part()

#from host import *

"""def modN(RP, savePathsInterval, systemClass, endTime
        , experimentName, opRunsFormula, mStepsPerOPRun
        , runsPerN, continueRun=False):

    startClock = time()
    maxWGSize = 512
    pathChanges = 0
    timestamp  = datetime.now().strftime("%Y-%m-%d/%H.%M.%S")
    filename = "results/" + experimentName + "/" + timestamp + "/data"

    # Create directory to store results and paths
    if not os.path.exists("results/" + experimentName + "/" + timestamp):
        os.makedirs("results/" + experimentName + "/" + timestamp)

    endN = RP.N
    startN = 8
    RP.S = 1

    if continueRun:
        startN = continueRun[1]
        RP.S = continueRun[2]

    # RP.N = startN, 16, .... , endN
    for RP.N in [2 ** i for i in range(int(np.log2(startN) + 0.5)
                 , int(np.log2(endN) + 1.0))]:

        # Make sure S is large enough to not use too many walkers in a WG
        while RP.nbrOfWalkersPerWorkGroup * RP.N / (2 ** RP.S) > maxWGSize:
            RP.S += 1


        # TODO: Should be able to specify
        RP.operatorRuns = opRunsFormula(RP.N, RP.S)
        RP.metroStepsPerOperatorRun = mStepsPerOPRun(RP.N, RP.S)
        RP.returnOperator = True

        KE = loadKernel(systemClass, RP)

        # If not first run create a new path by interpolating
        if not RP.N == startN:
            newPaths = np.zeros((RP.nbrOfWalkers, RP.N * systemClass.DOF))
            for walker in range(RP.nbrOfWalkers):
                for DOF in range(systemClass.DOF):
                    secondIndices = DOF * RP.N + np.array(range(0, RP.N / 2))
                    # Nodes from the old path is copied to every second node in
                    # new path.
                    newPaths[walker, secondIndices] = RKR.paths[walker
                            , secondIndices - RP.N * DOF / 2]

                    # Linear interpolation to create new nodes in between nodes
                    # from old path.
                    newPaths[walker, secondIndices + 1] = (newPaths[walker,
                            secondIndices] + newPaths[walker, secondIndices + 2]) / 2.0

        # If first run
        else:
            if not continueRun:
                # Create a random path
                initialPaths = (2.0 * np.random.rand(RP.nbrOfWalkers,RP.N *
                    systemClass.DOF) - 1.0) * 0.1
            else:
                initialPaths = np.zeros((RP.nbrOfWalkers, RP.N * systemClass.DOF))
                reader = csv.reader(open(continueRun[0]), delimiter='\t')
                k = 0
                for row in reader:
                    initialPaths[k] = [float(v) for v in row]
                    k += 1

        KE.paths.data.release()
        KE.paths = cl.array.to_device(KE.queue, initialPaths.astype(np.float32))

        RP.returnPaths = False
        nRuns = 1
        runsThisN = runsPerN(RP.N, RP.S)
        while (nRuns <= runsThisN or RP.N == endN):
            if time() - startClock > endTime:
                return

            # Save paths
            if nRuns % savePathsInterval == 0 or nRuns == runsThisN :
                print("Saving paths...")

                RP.returnPaths = True
                KE = loadKernel(systemClass, RP)
                RKR = runKernel(KE)
                nRuns += 1
                RP.returnPaths = False
                pathChanges += RP.getMetroStepsPerRun()
                output(filename, RP.N, time() - startClock, pathChanges
                        , RKR.acceptanceRate, RKR.operatorMean
                        , RP.beta, RP.S)
                f = open("results/" + experimentName + "/" + timestamp +
                         "/pathsN" + str(RP.N) + "episode" +
                         str(int(nRuns / savePathsInterval)), 'wb')
                csvWriter = csv.writer(f, delimiter='\t')
                for aPath in RKR.paths:
                    csvWriter.writerow(aPath)
                f.close()
                print("Paths saved!")

            RKR = runKernel(KE)
            nRuns += 1
            pathChanges += RP.getMetroStepsPerRun()
            output(filename, RP.N, time()-startClock, pathChanges
                   , RKR.acceptanceRate, RKR.operatorMean, RP.beta, RP.S)

        # Last run for this N so need to save paths
        if RP.N != endN:
            RP.returnPaths = True
            KE = loadKernel(systemClass, RP)
            RKR = runKernel(KE)
            RP.returnPaths = False
            pathChanges += RP.getMetroStepsPerRun()
            output(filename, RP.N, time()-startClock, pathChanges
                   , RKR.acceptanceRate, RKR.operatorMean, RP.beta, RP.S)

        # Change S to move acceptance rate in the right direction
        if RKR.acceptanceRate > 0.5:
            RP.S = min(i - 1, RP.S + 2)
        if 0.2 < RKR.acceptanceRate < 0.5:
            RP.S = min(i - 1, RP.S + 1)
        if RKR.acceptanceRate < 0.1:
            RP.S = max(1, RP.S - 1)


#def output(filename, N, t, pathChanges, acceptanceRate
           , operatorMean, beta, S):

    print("N: " + str(N) + "\tS: " + str(S) + "\tbeta: " +
          str(beta)+"\tAR: " + str(acceptanceRate) + "\tOP:s" +
          str(np.mean(operatorMean, axis = 1)))

    f_data = open(filename + ".tsv", 'a')

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
    f_data.close()
"""

def plotit(systemClass, dataFile, burnRatio, nOperators):

    # Physical constants
    a0 = 0.52917721092e-10  # Bohr radius in meters
    kB = 1.3806488e-23  # Boltzmanns constant in J/K

    # Conversion factors from units used in physical system to more
    # standard units.
    en2Kelv = sys.potentialUnit / kB
    len2BRad = sys.xUnit / a0

    # Load the data file and extract the parameters and op values.
    data = np.loadtxt(dataFile)
    generalData = data[:, 0:6]
    opData = data[:, 6:]

    nWalkers = len(opData[0]) / nOperators
    nRuns = len(generalData)

    N = data[:, 0].astype('int')
    S = data[:, 5].astype('int')
    AR = data[:, 3].astype('int')

    # Format the opData in a 3 dimensional array where the first index is the
    # operator, second is the walker, and the third is the run
    operatorMeansPerRun = np.empty((nOperators, nWalkers, nRuns))
    for i in range(nOperators):
        operatorMeansPerRun[i] = (np.array(opData[:, i * nWalkers:
                                 (i + 1) * nWalkers])).transpose()

    n = N[0]
    end = 0
    operatorMeansPerRunPerN = []
    burntOperatorMeansPerRunPerN = []
    operatorMeansPerN = []
    operatorMeansStdPerN = []
    runIndexPerN = []
    while n <= N[-1]:
        start = end
        if n == N[-1]:
            end = len(N)
        else:
            while N[end] == N[start]:
                end = end + 1
        runIndexPerN.append(range(start, end))
        operatorMeansPerRunPerN.append(operatorMeansPerRun[:, :, start:end])
        burntOperatorMeansPerRunPerN.append(operatorMeansPerRun[:, :, (start +
                                   int(float(end - start) * burnRatio)):end])
        operatorMeansPerN.append(burntOperatorMeansPerRunPerN[-1].mean(axis=2))
        n *= 2

    ####################
    # PLOT MEAN RADIUS
    ####################

    plt.figure(1)
    for i in range(nWalkers):
        plt.plot(operatorMeansPerRun[2, i, :] * sys.xUnit / a0)

    plt.xlabel('Run number')
    plt.ylabel('Mean radius [Bohr radius]')
    plt.grid(True)

    ####################
    # PLOT MEAN ENERGY
    ####################

    plt.figure(2)
    for i in range(len(operatorMeansPerRunPerN)):
        plt.plot(runIndexPerN[i],
                operatorMeansPerRunPerN[i][0, :, :].mean(axis=0) * en2Kelv,
                label=str(N[runIndexPerN[i][0]]))
    plt.xlabel('Run number')
    plt.ylabel('Energy, [K]')
    plt.legend(loc=4)
    plt.grid(True)

    ####################
    # PRINT CI FOR OPS
    ####################

    print('N= ' + str(N[-1]) + ", 95% CI:")

    mean = operatorMeansPerN[-1][0, :].mean() * en2Kelv
    std = operatorMeansPerN[-1][0, :].std() / np.sqrt(nWalkers) * en2Kelv
    print("Mean energy: [" + "%.4f" % (mean - 1.96 * std)
            + ',' + "%.4f" % (mean + 1.96 * std) + '] K')

    derivative = 0.5 / np.sqrt(operatorMeansPerN[-1][1, :].mean())
    mean = np.sqrt(operatorMeansPerN[-1][1, :].mean()) * len2BRad
    std = (operatorMeansPerN[-1][1, :].std() / np.sqrt(nWalkers)
            * derivative * len2BRad)
    print("Root mean square radius: [" + "%.4f" % (mean - 1.96 * std)
            + ',' + "%.4f" % (mean + 1.96 * std) + '] a0')

    mean = operatorMeansPerN[-1][2, :].mean() * len2BRad
    std = operatorMeansPerN[-1][2, :].std() / np.sqrt(nWalkers) * len2BRad
    print("Mean radius: [" + "%.4f" % (mean - 1.96 * std)
            + ',' + "%.4f" % (mean + 1.96 * std) + '] a0')

    plt.show()


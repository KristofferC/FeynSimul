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


import matplotlib.pyplot as plt
import csv
import numpy as np

def plotIt(systemClass, dataFile, burnRatio, nOperators):

    # Physical constants
    a0 = 0.52917721092e-10  # Bohr radius in meters
    kB = 1.3806488e-23  # Boltzmanns constant in J/K

    # Conversion factors from units used in physical system to more
    # standard units.
    en2Kelv = systemClass.potentialUnit / kB
    len2BRad = systemClass.xUnit / a0

    # Load the data file and extract the parameters and op values.
    data = np.loadtxt(dataFile)
    generalData = data[:, 0:5]
    opData = data[:, 5:]

    nWalkers = len(opData[0]) / nOperators
    nRuns = len(generalData)

    N = data[:, 0].astype('int')
    S = data[:, 4].astype('int')
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
        plt.plot(operatorMeansPerRun[2, i, :] * systemClass.xUnit / a0)

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

if __name__ == '__main__':
    print("This can no longer be run as a script.")

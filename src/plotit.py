import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import time


#plotit(systemClass, dataFile)
############################what to present##################
filename='5/modNBisect.tsv'
nbrOfWalkers = 64
nbrOfOperators = 3

sys.path.append(sys.path[0] + "/../../src/physical_systems/")
from lm2m2_3part import *
sys=Lm2m2_3part()


burnRatio=0.3

# Physical constants
a0=0.52917721092e-10  # Bohr radius in meters
kB = 1.3806488e-23 # Boltzmanns constant in J/K

# Conversion factors from units used in physical system to more standard units.
en2Kelv = sys.potentialUnit / kB
len2BRad = sys.xUnit / a0

# Load the data file and extract the parameters and data
data = np.loadtxt(filename)
generalData = data[:, 0:6]
opData = data[:, 6:]

# FIXA
nbrOfRuns = len(generalData)
if len(opData[0]) != nbrOfWalkers * nbrOfOperators:
    raise Exception('Wrong number of walkers or operators')

N = data[:, 0].astype('int')
S = data[:, 5].astype('int')
AR = data[:, 3].astype('int')

# Format the opData in a 3 dimensional array where the first index is the
# operator, second is the walker, and the third is the run
operatorMeansPerRun = np.empty((nbrOfOperators, nbrOfWalkers, nbrOfRuns))
for i in range(nbrOfOperators):
    operatorMeansPerRun[i] = (np.array(opData[:, i * nbrOfWalkers:(i + 1)
                                * nbrOfWalkers])).transpose()


n=N[0]
end=0
operatorMeansPerRunPerN = []
burntOperatorMeansPerRunPerN = []
operatorMeansPerN = []
operatorMeansStdPerN = []
runIndexPerN = []
while n<=N[-1]:
    start=end
    if n==N[-1]:
        end=len(N)
    else:
        while N[end]==N[start]:
            end=end+1
    runIndexPerN.append(range(start,end))
    operatorMeansPerRunPerN.append(operatorMeansPerRun[:,:,start:end])
    burntOperatorMeansPerRunPerN.append(operatorMeansPerRun[:,:,(start+int(float(end-start)*burnRatio)):end])
    operatorMeansPerN.append(burntOperatorMeansPerRunPerN[-1].mean(axis=2))
    n*=2


####################
# PLOT MEAN RADIUS
####################

plt.figure(1)
for i in range(nbrOfWalkers):
    plt.plot(operatorMeansPerRun[2, i, :] * sys.xUnit / a0)

plt.xlabel('Run number')
plt.ylabel('Mean radius [Bohr radius]')
plt.grid(True)


####################
# PLOT MEAN ENERGY
####################

plt.figure(2)
for i in range(len(operatorMeansPerRunPerN)):
    plt.plot(runIndexPerN[i],operatorMeansPerRunPerN[i][0,:,:].mean(axis=0) *
            en2Kelv, label=str(N[runIndexPerN[i][0]]))
plt.xlabel('Run number')
plt.ylabel('Energy, [K]')
plt.legend(loc = 4)
plt.grid(True)


print('N: ' + str(N[-1]) + , " 95% CI.")

mean = operatorMeansPerN[-1][0,:].mean() * en2Kelv
std = operatorMeansPerN[-1][0,:].std() / np.sqrt(nbrOfWalkers) * en2Kelv
print("Mean energy: [" + "%.4f" % (mean - 1.96 * std)
        + ',' + "%.4f" % (mean + 1.96 * std) + '] K')

derivative = 0.5 / np.sqrt(operatorMeansPerN[-1][1,:].mean())
mean=np.sqrt(operatorMeansPerN[-1-i][1,:].mean()) * len2BRad
std=operatorMeansPerN[-1-i][1,:].std()/np.sqrt(nbrOfWalkers) * derivative * len2BRad
print("Root mean square radius: ["+str(mean-2.0*std)+','+str(mean+2.0*std)+'] a0')

mean=operatorMeansPerN[-1][2, :].mean() * len2BRad
std=operatorMeansPerN[-1][2, :].std()/np.sqrt(nbrOfWalkers) * len2BRad
print("Mean radius: ["+str(mean - 1.96 * std)+','+str(mean + 1.96 * std)+'] a0')

plt.show()

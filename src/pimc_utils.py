# Does a run with the lm2m2 potential for 3 particles.

import sys
import numpy as np
from datetime import datetime
from time import time
import csv





# Time to run simul

# How often to save paths.


def modN(runparams, endTime, savePathsInterval, systemClass, endTime, endN,
        continueRun=False):

    startClock=time()
    maxWgSize = 512
    pathChanges = 0
    timestamp  = datetime.now().strftime("%Y-%m-%d/%H.%M.%S")
    filename = "results/"+experimentName+"/" + timestamp + "/modNBisect"

    # Create directory to store results and paths
    if not os.path.exists("results/"+experimentName+"/" + timestamp):
        os.makedirs("results/"+experimentName+"/" + timestamp)

    startN = 8
    RP.S = 1

    if continueRun:
        startN = continueRun[1]
        RP.S = continueRun[2]

    # RP.N = startN, 16, .... , endN
    for RP.N in [2 ** i for i in range(int(np.log2(startN) + 0.5))
                 , int(np.log2(endN) + 0.5))]:

        if time() - startClock > endTime:
            break

        # Make sure S is large enough to not use too many walkers in a WG
        while RP.nbrOfWalkersPerWorkGroup * RP.N / (2 ** RP.S) > maxWGSize:
            RP.S+=1

        RP.returnPaths = True
        RP.operatorRuns = max(RP.N / 8, 1)
        RP.metroStepsPerOperatorRun = 10
        RP.returnOperator = True

        KE = loadKernel(systemClass, RP)

        # If not first run create a new path by interpolating
        if RP.N != startN:
            newPaths = np.zeros((RP.nbrOfWalkers, RP.N * systemClass.DOF))
            for walker in range(RP.nbrOfWalkers):
                for DOF in range(systemClass.DOF):
                    secondIndices = DOF * RP.N + np.array(range(0, RP.N / 2))
                    # Nodes from the old path is copied to every second node in
                    # new path.
                    newPaths[walker, secondIndices]
                          = RKR.paths[walker , secondIndices - RP.N * DOF / 2]

                    # Linear interpolation to create new nodes in between nodes
                    # from old path.
                    newPaths[walker, secondIndices + 1] = (newPaths[threadID,
                            secondIndices] + newPaths[walker, secondIndices + 2]) / 2.0

            KE.paths.data.release()
        # If first run
        else:
            if not continueRun:
                # Create a random path
                initialPaths = (2.0 * np.random.rand((RP.nbrOfWalkers,RP.N *
                    systemClass.DOF)) - 1.0) * 0.1
            else:
                initialPaths = np.zeros((RP.nbrOfWalkers, RP.N * systemClass.DOF))
                reader = csv.reader(open(continueRun[0]), delimiter='\t')
                k = 0
                for row in reader:
                    initialPaths[k] = [float(v) for v in row]
                    k += 1

        KE.paths = cl.array.to_device(KE.queue,initialPaths.astype(np.float32))

        RP.returnPaths = False

        while (nRuns < max(RP.N / 8, 10) or RP.N == endN):
            if time() - startClock < endTime:
                return

            # Save paths
            if nRuns % savePathsInterval == 0:
                print("Saving paths...")

                RP.returnPaths = True
                RKR = runKernel(KE)
                RP.returnPaths = False
                pathChanges += RP.getMetroStepsPerRun()
                output(filename, RP.N, time() - startClock, pathChanges
                        , RKR.acceptanceRate,RKR.operatorMean,RP.beta,RP.S)
                f = open("results/"+experimentName+"/" + timestamp+"/pathsN"+str(RP.N)+"episode"+str(int(j/savePathsInterval)),'wb')
                csvWriter = csv.writer(f, delimiter='\t')
                for aPath in RKR.paths:
                    csvWriter.writerow(aPath)
                f.close()
                print("Paths saved!")

            RKR = runKernel(KE)
            pathChanges += RP.getMetroStepsPerRun()
            output(filename, RP.N,time()-startClock, pathChanges, RKR.acceptanceRate,RKR.operatorMean,RP.beta,RP.S)
            nRuns += 1

        # Last run for this N so need to save paths
        if RP.N != endN:
            RP.returnPaths = True
            RKR = runKernel(KE)
            pathChanges += RP.getMetroStepsPerRun()
            output(filename,RP.N,time()-startClock,pathChanges,RKR.acceptanceRate,RKR.operatorMean,RP.beta,RP.S)
T
        # Change S to move acceptance rate in the right direction
        if RKR.acceptanceRate > 0.5:
            RP.S = min(i - 1, RP.S + 2)
        if 0.2 < RKR.acceptanceRate < 0.5:
            RP.S = min(i - 1, RP.S + 1)
        if RKR.acceptanceRate < 0.1:
            RP.S = max(1, RP.S - 1)


def output(filename, N, t, pathChanges, acceptanceRate
           , operatorMean, beta, S, firstRun=False):

    print("N: " + str(N)+"\tS: "+str(S)+"\tbeta: "+str(beta)+"\tAR: " + str(acceptanceRate))

    f_data = open(filename + ".tsv",'a')


    if firstRun: # Add a heading describing the columns.
        f_data.write("#N\tTime\tpathChanges\tAR\tS")
         for i in range(len(operatorMean)):
             f_data.write("\tOperator " + str(i))

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

# modN("directory to last path").

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

sys.path.append(sys.path[0] + "/../src/physical_systems/")
from pisa3petter import *
from pisa4 import *
from pisa import *
from lm2m2_3part import *

experimentName="lm2m2_3part/"

RP = RunParams()
RP.nbrOfWalkers = 64
RP.N = 1024 * 64
RP.alpha = 0.5
RP.getOperator = True
RP.enablePathShift = False
RP.enableBisection = True
RP.enableSingleNodeMove = False
RP.enableGlobalPath = True
RP.enableGlobalOldPath = True
RP.enableParallelizePath = True
RP.returnBinCounts = False
RP.beta = 2000
RP.nbrOfWalkersPerWorkGroup=4

saveOpValNum = 0

# Time to run simul
endTime = 60*60*24*14

# How often to save paths.
savePathsInterval = 3000
systemClass = Lm2m2_3part()

RP.operators=(systemClass.energyOp,systemClass.meanSquaredRadiusOp,systemClass.meanRadiusOp)

def output(filename,N,t,pathChanges,acceptanceRate,operatorMean,beta,S):
    print("N="+str(N)+"\tS="+str(S)+"\tbeta="+str(beta)+"\tAR=" + str(acceptanceRate))
    #print("En="+str(operatorMean[0]*systemClass.potentialUnit/1.3806503e-23)+"("+str(systemClass.groundStateEnergy)+")")
    #print("Root mean sqr radius="+str(math.sqrt(operatorMean[1])*1e10*systemClass.xUnit)+"("+str(math.sqrt(systemClass.meanSquaredRadius))+")")
    #print("Mean rad="+str(operatorMean[2]*1e10*systemClass.xUnit)+"("+str(math.sqrt(systemClass.meanSquaredAtomDist))+")")
    f_data = open(filename + ".tsv",'a')
    f_data.write(str(N)+"\t")
    f_data.write(str(t)+"\t")
    f_data.write(str(pathChanges)+"\t")
    f_data.write(str(acceptanceRate)+"\t")
    f_data.write(str(beta)+"\t")
    f_data.write(str(S)+"\t")
    for walkerOperators in operatorMean:
        for op in walkerOperators:
            f_data.write(str(op)+"\t")

    #for tid in range(0,min(len(operatorValues),saveOpValNum)):
    #    f_data.write(str(operatorValues[tid])+"\t")
    f_data.write("\n")
    f_data.close()

def modN(continueRun=False):
    timestamp  = datetime.now().strftime("%Y-%m-%d/%H.%M.%S")
    if not os.path.exists("results/"+experimentName+"/" + timestamp):
        os.makedirs("results/"+experimentName+"/" + timestamp)

    startClock=time()
    pathChanges=0
    filename = "results/"+experimentName+"/" + timestamp + "/modNBisect"
    maxWGSize=1024
    endI=int(np.log2(RP.N)+0.5)
    RP.operatorRuns = 1
    RP.metroStepsPerOperatorRun=0
    RP.returnPaths=True;
    RP.returnOperator=True;
    if not continueRun:
        startI=3
        RP.S=1
        RP.N=2**startI
        while RP.nbrOfWalkersPerWorkGroup*RP.N/(2**RP.S)>maxWGSize:
            RP.S+=1
        i=startI
    else:
        i=continueRun[1]
        RP.S=continueRun[2]
        startI=continueRun[1]
    while i<endI+1 and time()-startClock<endTime:
        sleep(2)
        RP.N=2**i
        while RP.nbrOfWalkersPerWorkGroup*RP.N/(2**RP.S)>maxWGSize:
            RP.S+=1
        RP.operatorRuns = max(2**10/2**RP.S, 1)
        RP.metroStepsPerOperatorRun=10#max(16*64/((2**RP.S)),1)
        KE = loadKernel(systemClass, RP)

        if i>startI :
            newPaths=np.zeros((RP.nbrOfWalkers, RP.N * systemClass.DOF))
            for tid in range(0,RP.nbrOfWalkers):
                for particle in range(0,systemClass.DOF):
                    newPaths[tid,particle*RP.N+np.array(range(0,RP.N,2))]= RKR.paths[tid,particle*RP.N/2+np.array(range(0,RP.N/2))]
                    newPaths[tid,particle*RP.N+np.array(range(0,RP.N,2))+1]= \
                        (RKR.paths[tid,particle*RP.N/2+np.array(range(0,RP.N/2))]+RKR.paths[tid,particle*RP.N/2+(np.array(range(0,RP.N/2))+1)%(RP.N/2)])/2.0
            KE.paths.data.release()
            KE.paths = cl.array.to_device(KE.queue, newPaths.astype(np.float32)  )
        else:
            if not continueRun:
                initialPaths=np.zeros((RP.nbrOfWalkers,RP.N * systemClass.DOF))
                for wid in range(0,RP.nbrOfWalkers):
                    for particle in range(0,systemClass.DOF):
                        initialPaths[wid,particle*RP.N:(particle+1)*RP.N]=np.ones(RP.N)*(np.random.rand(1)*2.0-1.0)*0.1
                    KE.paths = cl.array.to_device(KE.queue,initialPaths.astype(np.float32)  )
            else:
                initialPaths=np.zeros((RP.nbrOfWalkers,RP.N * systemClass.DOF))
                reader = csv.reader(open(continueRun[0]), delimiter='\t')
                k=0
                for row in reader:
                    initialPaths[k]=[float(v) for v in row]
                    k+=1
                KE.paths = cl.array.to_device(KE.queue,initialPaths.astype(np.float32)  )

        RP.returnPaths = False
        print(KE.getStats())
        j=0
        while (j < max(RP.N/8,10) or i==endI ) and time()-startClock<endTime:
            if int(j/savePathsInterval)!=int((j-1)/savePathsInterval):
                print("Saving paths...")
                RP.returnPaths=True
                RKR = runKernel(KE)
                pathChanges+=RP.getMetroStepsPerRun()
                output(filename,RP.N,time()-startClock,pathChanges,RKR.acceptanceRate,RKR.operatorMean,RP.beta,RP.S)
                RP.returnPaths=False
                f=open("results/"+experimentName+"/" + timestamp+"/pathsN"+str(RP.N)+"episode"+str(int(j/savePathsInterval)),'wb')
                csvWriter=csv.writer(f, delimiter='\t')
                for aPath in RKR.paths:
                    csvWriter.writerow(aPath)
                f.close()
                print("Paths saved...")
            RKR = runKernel(KE)
            pathChanges+=RP.getMetroStepsPerRun()
            output(filename,RP.N,time()-startClock,pathChanges,RKR.acceptanceRate,RKR.operatorMean,RP.beta,RP.S)
            j+=1

        if i!=endI and time()-startClock<endTime:
            RP.returnPaths=True
            RKR = runKernel(KE)
            pathChanges+=RP.getMetroStepsPerRun()
            output(filename,RP.N,time()-startClock,pathChanges,RKR.acceptanceRate,RKR.operatorMean,RP.beta,RP.S)
        i+=1
        if RKR.acceptanceRate>0.5:
            RP.S+=1
        if RKR.acceptanceRate>0.2:
            RP.S+=1
        if RKR.acceptanceRate<0.1:
            RP.S-=1
        if RP.S>=i:
            RP.S=i-1
        if RP.S<1:
            RP.S=1

# modN("directory to last path").
modN()

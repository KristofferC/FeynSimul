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
        simTime=-1, verbosity=0, cont_S=None, cont_N=None, cont_path=None):
    """
    A function to automate some of the things one would commonly like to do
    when doing a simulation. This function helps with saving data to a file,
    resuming a simulation that was interrupted, increasing convergence speed as
    well as helps determining what ``N`` is needed for the Trotter product
    formula to be be accurate.

    This is done by starting with a low ``N`` and running the simulation for
    this ``N`` until convergence. A new path with twice as many nodes is then
    created by interpolating between the nodes for the old path. This procedure
    is done until the values of the operators does not noticable change when
    ``N`` is increased.

    The data for every run is saved to a file which can then later be used for
    post processing. This post processing can be done while the
    simulation is still going on to check for convergence and error intervals.

    The paths are saved at specified intervals and if the simulation crash it
    can be resumed from one of these saved path files. The paths are also
    automatically saved when ``N`` is increased.

    This function only works with the bisection sampling method right now.

    :type ka: :class:`kernel_args.kernelArgs` class
    :param ka: An instance of kernelArgs describing what kind of
               kernel to build.
    
    :type startXList: ndarray
    :param startXList: The initial paths for an individual simulation are all
                       straight lines. Where in space these straight lines will
                       lie are determined by the numbers in this array. One
                       number for each path.
                       These could for example be uniform random numbers to get
                       a nice spread of initial paths between different
                       independent simulations.

    
    :type savePathsInterval: int
    :param savePathsInterval: How many runs to wait before saving the path to
                              disk.

    :type experimentName: string
    :param experimentName: A label for the simulation. The results will be
                           stored in a map with this name.

    :type opRunsFormula: function
    :param opRunsFormula: f(N, S) that returns the number of operator runs the
                          kernel should do before returning. For an
                          approximately constant run time this should be
                          inversely proportional to :math:`2^S`.

    :type mStepsPerOPRun: function
    :param mStepsPerOPRun: f(N, S) that returns the number of metropolis
                           samples to do before calculating the operator.

    :type runsPerN: function
    :param runsPerN: f(N, S) that returns the number of runs to do before
                     increasing the number of nodes in the path.

    :type maxWGSize: int
    :param maxWGSize: Maximum work group size for the GPU for this specific
                      kernel... TODO: automate this with
                      cl.kernel_work_group_info......

    :type continueRun: boolean
    :param continueRun: Enable to indicate that this will continue a simulation
                        that was 
  
    :type cont_path: string
    :param cont_path: The path to a file where a path from another simulation
                      was saved.
    
    :type cont_S: int
    :param cont_S: The value of ``S`` to start using when continuing an earlier
                   run.

    :type cont_N: int
    :param cont_N: The value of ``N`` to start using when continuing an earlier
                   run.

    :type finalN: int
    :param finalN: At what point the number of nodes in the path should stop
                  increasing. Defaulted to -1 which means it will never stop.

    :type simTime: float
    :param simTime: The total time in seconds that the simulation will run for. 
                    Defaulted to -1 which means that it will never stop.

    :type verbosity: int
    :param verbosity: Prints out extra information to stdout if greater than 0.
    """

    kaMod = copy.copy(ka)

    if kaMod.enableBisection == False:
        raise Exception('This function currently only works with bisection'
                         ' sampling')

    if continueRun:
        if cont_path == None:
            raise NameError('cont_path need to be defined if continueRun'
                            'is enabled')
        if cont_S == None:
            raise NameError('cont_S need to be defined if continueRun'
                            ' is enabled')
        if cont_N == None:
            raise NameError('cont_N need to be defined if continueRun'
                            ' is enabled')

        kaMod.N = cont_N
        kaMod.S = cont_S
    else:
        kaMod.N = ka.N
        kaMod.S = 1
    startClock = time()
    pathChanges = 0
    timestamp  = datetime.now().strftime("%Y-%m-%d/%H.%M.%S")
    filename = "results/" + experimentName + "/" + timestamp + "/data"

    # Create directory to store results and paths
    if not os.path.exists("results/" + experimentName + "/" + timestamp):
        os.makedirs("results/" + experimentName + "/" + timestamp)
    
   

    
    firstNRun = True

    while finalN == -1 or kaMod.N <= finalN:
        # Make sure S is large enough to not use too many walkers in a WG
        #while ka.nbrOfWalkersPerWorkGroup * kaMod.N / (2 ** kaMod.S) > maxWGSize:
        #    kaMod.S += 1
        minS=int(np.ceil(np.log2(ka.nbrOfWalkersPerWorkGroup *
            kaMod.N/float(maxWGSize))))
        maxS=int(np.log2(kaMod.N))
        if maxS<minS:
            raise Exception('Too many walkers per workgroup')
        kaMod.S=min(maxS,max(minS,kaMod.S))
        kaMod.operatorRuns = opRunsFormula(kaMod.N, kaMod.S)
        kaMod.metroStepsPerOperatorRun = mStepsPerOPRun(kaMod.N, kaMod.S)

        kernel=PIMCKernel(kaMod, verbose = verbosity >= 2)
        if verbosity>0:
            print("New kernel compiled, getStats():")
            print(kernel.getStats())
            print("Run results:")

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
                reader = csv.reader(open(cont_path), delimiter='\t')
                k = 0
                for row in reader:
                    initialPaths[k] = map(float, row)
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
            if simTime != -1 and time() - startClock > simTime:
                if verbosity>0:
                    print("Time limit reached!")
                return
            kernel.run()
            pathChanges += kernel.getMetroStepsPerRun()
            _output(filename, kaMod.N, time()-startClock, pathChanges
                   , kernel.getAcceptanceRate(), kernel.getOperators()
                   , ka.beta, kaMod.S,verbose=verbosity>0)
            # Save paths
            if nRuns % savePathsInterval == 0 or nRuns == runsThisN :
                with open("results/" + experimentName + "/" + timestamp +
                         "/pathsN" + str(kaMod.N) + "episode" +
                         str(int(nRuns / savePathsInterval)), 'wb') as f:
                    csvWriter = csv.writer(f, delimiter='\t')
                    for aPath in kernel.getPaths():
                        csvWriter.writerow(aPath)
                if verbosity>0:
                    print("Paths saved!")
            nRuns += 1

        with open("results/" + experimentName + "/" + timestamp +
                 "/timestamps", 'ab') as f2:
            f2.write(str(kaMod.N) + ': ' + str(time()-startClock) + '\n')
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
            kaMod.N *= 2

def _output(filename, N, t, pathChanges, acceptanceRate
           , operatorMean, beta, S, verbose=False):
    """
    Used by the modN function to save data from a run into a file in ASCII
    format.

    :type: filename: string
    :param filename: path + filename of the data file.

    :type N: int
    :param N: Total nodes in a path for the current run.

    :type t: float
    :param t: Total time elapsed since start of simulation.

    :type pathChanges: int
    :param pathChanges: Total number of changes done to a path since the start
                        of simulation

    :type acceptanceRate: float
    :param acceptanceRate: The acceptance rate for this run.

    :type operatorMean: ndarray
    :param operatorMean: Mean of the operator calculated by each independent
                         simulations. Will be an array with the length as the
                         product of the number of operators and the number of
                         independent simulations

    :type beta: float
    :param beta: Euclidian time used for the simulation. Should be constant
                 between runs.

    :type S: int
    :param S: Parameter for bisection sampling.

    :type verbose: boolean
    :param verbose: If enabled prints a bit of info to stdout.
    """

    if verbose:
        print("N: " + str(N) + "\tS: " + str(S) + "\tbeta: " +
                str(beta)+"\tAR: " + str(acceptanceRate) + "\tOPs: " +
            str(np.mean(operatorMean, axis = 1)))

    with open(filename + ".tsv", 'a') as f_data:
        # Add a heading describing the columns if file is empty.
        if  os.path.getsize(filename + ".tsv") == 0:
            f_data.write("#Beta: " + str(beta) + "\n")
            f_data.write("#N\tTime\tpathChanges\tAR\tS")
            for i in range(len(operatorMean)):
                for j in range(len(operatorMean[0])):
                 f_data.write("\tOperator " + str(i) + ", Walker " + str(j) )
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

def plotModN(filename):
    """    
    Used to plot the result of result files generated by modN.

    :type: filename: string
    :param filename: path + filename of the data file.
    """

    reader = csv.reader(open(filename), delimiter='\t')
    cleanedData = [row for row in reader if not '#'==row[0][0]]
    
    N=np.array([int(row[0]) for row in cleanedData])
    AR=np.array([float(row[3]) for row in cleanedData])
    energy=np.array([float(row[6]) for row in cleanedData])
    meanSquareRadius=np.array([float(row[7]) for row in cleanedData])
    sys.path.append(sys.path[0] + "/../../FeynSimul/physical_systems/")

    plotData=meanSquareRadius
    n=N[0]

    end=0
    ni=0
    smoothedPlotData=[]
    while n<=N[-1]:
        start=end
        if n==N[-1]:
            end=len(plotData)
        else:
            while N[end]==N[start]:
                end=end+1
        smoothedPlotData.append(np.zeros((end-start,2)))
        smoothness=int(n/512)
        kernel=np.ndarray(smoothness*2+1)
        if smoothness==0:
            kernel[0]=1.0
        else:
            for i in range(smoothness*2+1):
                kernel[i]=np.exp(-((2.0*float(i-smoothness)/float(smoothness))**2.0))
        for i in range(start,end):
            tot=0.0
            for j in range(smoothness*2+1):
                if i+j-smoothness>start and i+j-smoothness<end:
                    smoothedPlotData[ni][i-start,0]=smoothedPlotData[ni][i-start,0]+kernel[j]*plotData[i+j-smoothness]
                    tot=tot+kernel[j]
            smoothedPlotData[ni][i-start,0]=smoothedPlotData[ni][i-start,0]/tot
            smoothedPlotData[ni][i-start,1]=i
        n=2*n
        ni=ni+1
    ##present data
    burnRatio=0.3
    useIndicies=[i for i in range(len(N)) if N[i]==N[-1]]
    endBurn=useIndicies[0]+int(float(useIndicies[-1]-useIndicies[0])*burnRatio)
    useIndicies=range(endBurn,useIndicies[-1]+1)
    """
    print("Meaned energy, last N-group: "+str(np.mean([energy[i] for i in useIndicies])*sys.potentialUnit/1.3806503e-23)+' K')
    a0=0.52917721092e-10
    print("Root mean square radius, last N-group:
    "+str(np.sqrt(np.mean([meanSquareRadius[i] for i in
    useIndicies]))*sys.xUnit/a0)+' a0')"""
    import matplotlib.pyplot as plt
    #plt.plot(plotData)
    for i in smoothedPlotData:
        plt.plot(i[:,1],i[:,0],label=str(N[i[0,1]]))
    plt.axis([0.0,len(plotData),-25.0,0.0])
    plt.grid(True)
    plt.legend(loc=4)
    plt.show()

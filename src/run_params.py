"""
This module defines a single class, C{RunParams} which is used to keep track
of parameters affecting a run.
"""
__docformat__ = 'epytext en' 


class RunParams:
    """
    Stores variables that are needed to perform a PIMC run.

    Variables that always need to be defined in the class before a run are:

      - B{N} (int) E{-} The number of nodes in a path.
      - B{beta} (float) E{-} Euclidian time difference between first and last node
      in the path.
      - B{nbrOfWalkers} (int) E{-} The number of paths that should be worked on
      simultaneously.
      - B{operatorRuns} (int) E{-} The number of times the operator should be
      calculated before returning the results from the GPU.
      - B{metroStepsPerOperatorRun} (int) E{-} The paths between operator
      measurements should be uncorrelated so it makes no sense to calculate
      determines how many changes should be made to the path before the operator
      should be calculated along the whole path.

    There are also options to make regarding what type of changes that should
    be made to the path. They are listed below together with their
    corresponding parameters.

      - B{enableBisection} (bool) E{-} Uses the bisection algorithm to make
      changes to the path. 
        - B{S} (int) E{-} The total bisection depth'
        - B{enableParallelizePath} (bool) E{-} This makes it so that many
        threads work on the same path simultaneously. The number of threads
        that will work on a path is B{N} / (2**B{S}).

      - B{enablePathShift} (bool) E{-} Takes an interval of the path and offsets the
      whole interval.
        - B{PSAlpha} (float) E{-} The offset of the interval of the path will
        be a random number between -PSAlpha and PSAlpha.

      - B{enableSingleNodeMove} (bool) E{-} Offsets a single node in the path.
        - B{alpha} (float) E{-} The offset of the node of the path will be a random
        number between E{-} B{alpha} and B{alpha}.

    Caclulation of the probability density is enabled with binsEnabled (bool).
    Variables that affect how the probability density is calculated are
      - B{xmin} (float) & B{xmax} (float) E{-} A hypercube with the sides B{xmin} 
      to B{xmax} where the probability density will be calculated inside.
      - B{binResolutionPerDOF} E{-} Determines how many discretization points to do
      in a degree of freedom. The total number of bins will be binResolutionPerDOF
      to the power of the degrees of freedom in the system.

    You can chose what type of memory the paths that the threads are working on
    should be stored with B{enableGlobalPath} (bool). With enableGlobalPath they are
    stored in the global memory and else in the local memory.

    If B{returnPaths} (bool) is enabled the paths are fetched from the GPU after
    a run.
    """
    def __init__(self, enableBins=False, returnBinCounts=False, enableOperator=True,enableCorrelator=False,
                 enablePathShift=False, enableBisection=False,
                 enableSingleNodeMove=False, enableParallelizePath=False,
                 enableGlobalPath=False, enableGlobalOldPath=False,
                 returnPaths=False):
        """
        Creates a new C{RunParams}. Sets most options to False.
        """

        self.enableBins = enableBins
        self.returnBinCounts = returnBinCounts
        self.enableOperator = enableOperator
        self.enableCorrelator = enableCorrelator
        self.enablePathShift = enablePathShift
        self.enableBisection = enableBisection
        self.enableSingleNodeMove = enableSingleNodeMove
        self.enableParallelizePath = enableParallelizePath
        self.enableGlobalPath = enableGlobalPath
        self.enableGlobalOldPath = enableGlobalOldPath
        self.returnPaths = returnPaths
        self.returnOperator=False
        self.returnCorrelator=False
        self.returnBinCounts=False

    def getMetroStepsPerRun(self):
        """
        @rtype: int
        @return: The total amount of metropolis steps that will be 
        done in a kernel run
        """
        if self.enableParallelizePath:
            nbrOfThreadsPerWalker = self.N / (2 ** self.S)
            return (self.nbrOfWalkers*nbrOfThreadsPerWalker*self.metroStepsPerOperatorRun * self.operatorRuns 
                *self.enableBisection)
        else:
            return (self.nbrOfWalkers*self.metroStepsPerOperatorRun * self.operatorRuns 
                * (self.enableBisection + self.enablePathShift	
                + self.enableSingleNodeMove))

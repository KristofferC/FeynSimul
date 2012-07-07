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

"""
This module defines a single class, :class:`KernelArgs` which is used to keep track
of parameters needed to create a Kernel.
"""

class KernelArgs:

    def __init__(self,
                 system,
                 nbrOfWalkers,
                 N,
                 beta,
                 enableOperator,
                 enableCorrelator,
                 enableBisection,
                 enablePathShift,
                 enableSingleNodeMove,
                 enableBins,
                 nbrOfWalkersPerWorkGroup,
                 enableDouble=None,
                 S=None,
                 alpha=None,
                 PSAlpha=None,
                 enableParallelizePath=None,
                 metroStepsPerOperatorRun=None,
                 enableGlobalPath = None,
                 enableGlobalOldPath = None,
                 operators=None,
                 operatorRuns=None,
                 correlators=None,
                 xMin=None,
                 xMax=None,
                 binResolutionPerDOF=None):
        """
        Creates an instance of a class that contains different parameters for the run.
        Those parameters that have default values does not necessarily need to be
        defined.

        :type system: :class:`physical_system` class instance
        :param system: The physical system to simulate.

        :type nbrOfWalkers: int
        :param nbrOfWalkers: The number of independent simulations that
                             are started simultaneously.

        :type N: int
        :param N: The number of nodes in a path.

        :type beta: float
        :param beta: Euclidian time difference between first and last node
                     in the path. 1/beta can be regarded as a measurement
                     of the temperature of the system.

        :type enableOperator: boolean
        :param enableOperator: Determines if the simulation should do
                               calculations on the operators defined.

        :type enableCorrelator: boolean
        :param enableCorrelator: Determines if the simulation should do
                                 calculations on the correlators defined.

        :type enableBisection: boolean
        :param enableBisection: Determines if the bisection sampling method
                                should be used.

        :type enablePathShift: boolean
        :param enablePathShift: Determince if the sampling method were a chunk
                                of the path is displaced in one sample step.

        :type enableSingleNodeMove: boolean
        :param enableSingleNodeMove: Determines if the sampling method were
                                     one node is displaced should be used.

        :type enableBins: boolean
        :param enableBins: Determines if the simulation should do calculations
                           on the probability density of the system.

        :type nbrOfWalkersPerWorkGroup: int
        :param nbrOfWalkersPerWorkGroup: Determines how many of the simulations 
                                         should be launched in one OpenCL 
                                         workgroup.

        :type S: int 
        :param S: Parameter to the bisection sampling algorithm. 2 ** S is the
                  number of nodes on sample step will update.

        :type alpha: float
        :param alpha: Parameter to the single node move sampling algorithm. A
                      node will be displaced between -alpha and alpha uniformly.

        :type PSAlpha: float
        :param PSAlpha: Parameter to the smapling algorithm that moves a chunk
                        of the path. The chunk will be displaced between
                        -PSAlpha and PSAlpha uniformly.

        :type enableParallelizePath: boolean
        :param enableParallelizePath: Parameter to the bisection sampling algorithm.
                                      Determines if many threads should work on the same
                                      simulations. This can greatly increase convergence
                                      speed.
        
        :type metroStepsPerOperatorRun: int
        :param metroStepsPerOperatorRun: The number of sampling iterations to do before
                                         the operator is calculated. Successive paths
                                         are correlated so calculating the operator
                                         too often is a waste.
       
        :type enableGlobalPath: boolean
        :param enableGlobalPath: Determines if the paths should be stored in local or
                                 global memory. Too large N and this might have to be
                                 enabled.
        
        :type enableGlobalOldPath: boolean
        :param enableGlobalOldPath: Determines if the path interval changed
                                    that might be reverted to if the Metropolis step is rejected
                                    will be stored in global memory.
        
        :type operators: tuple of strings
        :param operators: The operators to be calculated.
        
        :type operatorRuns: int
        :param operatorRuns: How many times the operator should be run before the kernel
                             returns.
        
        :type correlators: tuple of strings
        :param correlators: The correlators to be calculated
        
        :type xMin: float
        :param xMin: The "left" corner of the hypercube in which to calculate the
                     probability density.
        
        :type xMax: float
        :param xMax: The "right" corner of the hypercube in which to calculate the
                     probability density.
        
        :type binResolutionPerDOF: int
        :param binResolutionPerDOF: How finely the discretization of space is to be
                                    do`e for the probability density calculations.
                                    The total number of sub cubes are
                                    binResolutionPerDOF ** DOF.
                                    
        :type enableDouble: boolean
        :param enableDouble: Determine if double or single precision is to be
                             used for floats on the GPU. Default is single
                             precision.
	"""
        self.system = system
        self.nbrOfWalkers = nbrOfWalkers
        self.N = N
        self.S = S
        self.beta = beta
        self.nbrOfWalkersPerWorkGroup = nbrOfWalkersPerWorkGroup
        self.enablePathShift = enablePathShift

        self.enableSingleNodeMove = enableSingleNodeMove
       
        self.alpha = alpha

        self.enableOperator = enableOperator
     

        self.operatorRuns = operatorRuns
        self.operators = operators
        self.metroStepsPerOperatorRun = metroStepsPerOperatorRun

        self.enableCorrelator = enableCorrelator
      

        self.correlators = correlators

        self.enableBins = enableBins
       
        self.xMin = xMin
        self.xMax = xMax
        self.binResolutionPerDOF = binResolutionPerDOF


        self.enableBisection = enableBisection
        self.enableParallelizePath = enableParallelizePath
        self.enableGlobalPath = enableGlobalPath
        self.enableGlobalOldPath = enableGlobalOldPath
        self.enableDouble = enableDouble


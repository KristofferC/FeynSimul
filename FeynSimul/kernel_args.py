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
This module defines a single class, C{KernelArgs} which is used to keep track
of parameters needed to create a Kernel.
"""

class KernelArgs:

    def __init__(self,
                 system,
                 nbrOfWalkers,
                 N,
                 S,
                 beta,
                 enableOperator,
                 enableCorrelator,
                 enableBisection,
                 enablePathShift,
                 enableSingleNodeMove,
                 enableBins,
                 nbrOfWalkersPerWorkGroup,
                 alpha=None,
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


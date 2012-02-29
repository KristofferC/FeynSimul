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
        if enableSingleNodeMove:
            if alpha == None:
                raise NameError('alpha need to be set if enableSingleNodeMove'
                                ' is True.')
            self.alpha = alpha

        self.enableOperator = enableOperator
        if enableOperator:
            if operatorRuns == None:
                raise NameError('operatorRuns need to be set if enableOperator'
                                ' is True.')
            if operators == None:
                raise NameError('operators need to be set if enableOperator'
                                ' is True.')
            if metroStepsPerOperatorRun == None:
                raise NameError('metroStepsPerOperatorRun need to be set'
                                ' if enableOperator is True.')

            if not isinstance(operators, tuple):
                raise TypeError('operators need to be a tuple of strings.')

            self.operatorRuns = operatorRuns
            self.operators = operators
            self.metroStepsPerOperatorRun = metroStepsPerOperatorRun

        self.enableCorrelator = enableCorrelator
        if enableCorrelator:
            if correlators == None:
                raise NameError('correlators need to be set if enableCorrelator'
                                ' is True.')
            if not isinstance(correlators, tuple):
                raise TypeError('correlators need to be a tuple of strings.')

            self.correlators = correlators

        self.enableBins = enableBins
        if enableBins:
            if xMin == None:
                raise NameError('xMin need to be set if enableBins'
                                ' is True.')
            if xMax == None:
                raise NameError('xMax need to be set if enableBins'
                                ' is True.')
            if binResolutionPerDOF == None:
                raise NameError('binResolutionPerDOF need to be set if'
                                ' enableBins is True.')
            self.xMin = xMin
            self.xMax = xMax
            self.binResolutionPerDOF = binResolutionPerDOF


        self.enableBisection = enableBisection
        if enableBisection:
            if enableGlobalPath == None:
                raise NameError('enableGlobalPath need to be set if'
                                ' enableBisection is True.')
            if enableGlobalOldPath == None:
                raise NameError('enableGlobalOldPath need to be set'
                                 'if enableBisection is True.')
            if enableParallelizePath == None:
                raise NameError('enableParallelizePath need to be set'
                                'if enableBisection is True')
            self.enableParallelizePath = enableParallelizePath
            self.enableGlobalPath = enableGlobalPath
            self.enableGlobalOldPath = enableGlobalOldPath


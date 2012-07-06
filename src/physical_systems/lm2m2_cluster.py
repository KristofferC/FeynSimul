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

# -*- coding: utf-8 -*-

import jacobi
import partSys
import sympy
from sympy.utilities.codegen import codegen
    
class Lm2m2_cluster:
    def __init__(self,n,useJacobi=True,simplify=True):
        #kernel beta of 1000 gives real temp of 15 mK
        if useJacobi:
            self.DOF = (n-1)*3
        else:
            self.DOF = n*3
            raise Exception("disablement of jacobi isn't implemented yet")
        if n<2:
            raise Exception("n must be >=2")
        if 3<=n<=10:
            groundStates=[0.0872,0.387,0.901,1.605,2.478,3.489,4.641,5.904]#K
            self.groundStateEnergy=groundStates[n-3]
        self.xUnit= 1e-10 #Meter (xUnit is one Angstom)
        self.potentialUnit= 1.3806503e-23 #Joule (potentialUnit is one K)
        if simplify:
            r=[sympy.symbols('r'+str(i)) for i in range(3*n)] #cartesian coords
            rm=[sympy.Matrix([r[i*3],r[i*3+1],r[i*3+2]]) for i in range(n)]
                #real masses
            x=[sympy.symbols('x'+str(i)) for i in range(3*n)] #jacobi coords
            xm=[sympy.Matrix([x[i*3],x[i*3+1],x[i*3+2]]) for i in range(n)]
                #jacobi masses
            m=[sympy.symbols('m')]*n
            lm2m2sqrt=sympy.symbols('lm2m2sqrt')
            #self.symPotential = jacobi.toJacobiCoord([
            #    partSys.getHyperRadialPotential(rm,lambda
            #            a,b:((a-b).transpose()*(a-b))[0],[lm2m2sqrt])],rm,m,xm)[0].simplify()
            print('hello!!')
            print(partSys.getHyperRadialPotential)
            print('godbye!!')
            
            self.symEnergyOp = sum([sympy.diff(self.symPotential,xi)*xi/2 for xi in x])
            def toCode(expr,funName):
                [(sourceFilename,source),(headerFilenam,header)]=codegen((funName,expr),'C','codeGened',header=False)
                return source
            arglist=''
            for i in range(3*n):
                arglist+='x'+str(i)+(',' if i!=3*n else '')
            self.potential="potential("+arglist+")"
        else:
            raise Exception("disablement of simplify isn't implemented yet")
        self.userCode = toCode(sympy.symbols('x')**2,'potential')

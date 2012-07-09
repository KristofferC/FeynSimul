import jacobi
import partSys
import sympy
from sympy.utilities.codegen import codegen

class Lm2m2_cluster:
    def __init__(self,n,useJacobi=True,simplify=True):
        #kernel beta of 1000 gives real temp of 15 mK
        d=3 #we are living in 3 space-like dimensions
        if useJacobi:
            self.DOF = (n-1)*d
        else:
            self.DOF = n*d
            raise Exception("disablement of jacobi isn't implemented yet")
        if n<2:
            raise Exception("n must be >=2")
        if 3<=n<=10:
            groundStates=[0.0872,0.387,0.901,1.605,2.478,3.489,4.641,5.904]#K, ref?
            self.groundStateEnergy=groundStates[n-3]
        self.xUnit= 1e-10 #Meter (xUnit is one Angstom)
        self.potentialUnit= 1.3806503e-23 #Joule (potentialUnit is one K)
        if simplify:
            rm=[sympy.symbols(['r'+str(i*d+j) for j in range(d)]) for i in range(n)] #normal coords
            x=sympy.symbols(['x'+str(i) for i in range(n*d)])
            xm=[[x[i*d+j] for j in range(d)] for i in range(n)] #jacobi coords
            mm=[sympy.symbols('m')]*n
            dlm2m2sqrtdr=sympy.Symbol('dlm2m2sqrtdr')
            class lm2m2sqrt(sympy.core.function.Function):
                def fdiff(self, argindex=1):
                    return dlm2m2sqrtdr(self.args[0])
            self.symPotential = sympy.simplify(jacobi.toJacobiCoord([
                sympy.simplify(partSys.getHyperRadialPotential(rm,lambda
                    a,b:sum((a[i]-b[i])**2 for i in
                        range(d)),[lm2m2sqrt]))],rm,mm,xm)[0])
            self.symEnergyOp = self.symPotential+sympy.simplify(
                    sum(sympy.diff(self.symPotential,xi)*xi/2 for xi in x))
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

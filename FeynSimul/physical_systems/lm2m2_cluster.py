import jacobi
import partSys
import sympy
from sympy.utilities.codegen import codegen

class Lm2m2_cluster:
    def __init__(self,n,T,useJacobi=True,simplify=True,verbose=False):
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
        
        HeMass = 6.64647617e-27 #Kg, mass of one He atom
        hbar=1.05457148e-34 #m**2kg/s
        Kb=1.3806503e-23 #J/K 

        self.rUnit= 1e-10 #(before jacobi) Meter (xUnit is one Angstom)
        self.potentialUnit= Kb #Joule (potentialUnit is 1 K*Kb)
        self.beta=self.potentialUnit/(Kb*T) #This beta gives wanted T
        m=self.m=HeMass*self.rUnit**2/hbar**2*self.potentialUnit
        print m
        if simplify:
            rm=[sympy.symbols(['r'+str(i*d+j) for j in range(d)]) for i in
                    range(n)] #normal coords, in Angstrom
            x=sympy.symbols(['x'+str(i+1) for i in range(n*d)]) #we should really have indices start at 0 here...
            xm=[[x[i*d+j] for j in range(d)] for i in range(n)] #jacobi coords
            mm=[m]*n#[sympy.symbols('m')]*n
            dlm2m2sqrtdr=sympy.Symbol('dlm2m2sqrtdr')

            class lm2m2sqrt(sympy.Function):
                def fdiff(self, argindex=1):
                    return dlm2m2sqrtdr(self.args[0])
            
            def pow2ToSqr(s):  #doesn't work with nested () but we shouldn't
                i=-1           #have those
                while True:
                    i=s.find('pow(',i+1)
                    if i is -1: break
                    else:
                        k=s.find(')',i)
                        if s[k-3:k]==', 2':
                            s=s[:i]+'sqr'+s[i+3:k-3]+s[k:]
                return s
            def toCode(expr):
                [(sourceFilename,source),(headerFilenam,header)]=codegen(
                        ('funName',expr),'C','codeGened',header=False)
                code=pow2ToSqr(source).replace('double','float')
                code=code[code.find('return')+6:]
                code=code[:code.find(';')]
                return code
            
            if verbose: print('Making potential in space coordinates...')
            p=sympy.simplify(partSys.getHyperRadialPotential(rm,lambda
                                    a,b:sum((a[i]-b[i])*(a[i]-b[i]) for i in
                                                                range(d)),[lm2m2sqrt]))
            
            if verbose: print('Transforming to Jacobi coordinates...')
            self.symPotential = jacobi.toJacobiCoord([p],rm,mm,xm)[0]
            
            if verbose: print('Simplifying potential...')
            self.symPotential = sympy.simplify(self.symPotential)
            
            if verbose: print('Making virial energy estimator...')
            self.symEnergyOp = self.symPotential+(
                    sum(sympy.diff(self.symPotential,xi)*xi/2 for xi in x))
            
            if verbose: print('Making MSR operator...')
            self.symMSROp = sum((rm[i][di]-rm[j][di])**2 for di in range(3) for i in range(n) for j in range(i))
            
            if verbose: print('Transforming to jacobi coordinates...')
            self.symMSROp = sympy.simplify(jacobi.toJacobiCoord([self.symMSROp],rm,mm,xm)[0])
            
            if verbose: print('Generating C-code for potential...')
            self.potential=toCode(self.symPotential) 
            
            if verbose: print('Generating C-code for energy operator...')
            self.energyOp=toCode(self.symEnergyOp)

            if verbose: print('Generating C-code for MSR operator...')
            self.meanSquaredRadiusOp=toCode(self.symMSROp)
        else:
            raise Exception("disablement of simplify isn't implemented yet")
        with open('FeynSimul/physical_systems/lm2m2.h','r') as f:
            self.userCode = f.read()

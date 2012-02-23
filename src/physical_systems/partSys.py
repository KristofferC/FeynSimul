import jacobi

# x=particle positions, interactions is a list of two body interactions, three
# body interactions, ... . All expressed in relevant squared hyperradii.
def getHyperRadialPotential(x,distFun,interactions):
    _sum = lambda a: reduce(lambda b,c:b+c,a) #stupid python sum adds a zero
    def hyperSqrRadius(y):
        return _sum(distFun(y[i],y[j]) for i in range(len(y)) for j in range(i+1,len(y)))
    def act(l,interaction,part=[]):
        if len(part)==l:
            return interaction(hyperSqrRadius([x[p] for p in part]))
        return _sum(act(l,interaction,part+[j]) for j in
            range(len(x)-l+len(part)+1) if part==[] or j>part[-1])
    return _sum(act(i+2,interactions[i]) for i in range(len(interactions)))

''':
#just some testing of this code. It should be able to take any objects as x,
#interactions must return objects that together can be operands to +.
import sympy
n=4;d=2
x=[sympy.Matrix([sympy.symbols('x'+str(i)), sympy.symbols('y'+str(i))]) for i in range(n)]
m=[sympy.symbols('m')]*n
f=[(lambda shr,i=i:sympy.symbols('f'+str(i))(shr)) for i in range(n-1)]
print(f)
sympy.pprint( getHyperRadialPotential(x,lambda
    a,b:((a-b).transpose()*(a-b))[0],f).simplify()  )'''

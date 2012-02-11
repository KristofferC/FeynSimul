import jacobi

# x=particle positions, interactions is a list of two body interactions, three
# body interactions, ... . All expressed in relevant squared hyperradiuses.
def getPotential(x,interactions):
    def hyperSqrRadius(y):
        return sum([sum([(y[a]-y[b])**2 for b in range(len(y)) if a<b]) for 
            a in range(len(y))])
    def act(l,interaction,part=[]):
        if len(part)==l:
            return interaction(hyperSqrRadius([x[p] for p in part]))
        return sum([act(l,interaction,part+[j]) for j in range(len(x)) 
            if part==[] or j>part[-1]])
    return sum([act(i+2,interactions[i]) for i in range(len(interactions))])

import sympy
n=7
x=[sympy.symbols('x'+str(i)) for i in range(n)]
m=[sympy.symbols('m')]*n
f=sympy.symbols('f')
print(jacobi.toJacobiCoord([getPotential(x,
    [lambda sr:f(sr)])],x,m)[0].simplify())

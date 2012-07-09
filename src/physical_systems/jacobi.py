import sympy as sp

def jacobiCoord(x,m):
    jx=range(len(m))
    _sum = lambda a: reduce(lambda b,c:b+c,a)
    sm=[_sum(m[:i+1]) for i in range(len(m))]
    for i in range(len(m)):
        if i == 0:
            jm = m[0] * m[1] / (m[0] + m[1])
            jx[i] = (x[0] - x[1]) * sp.sqrt(jm)
        elif i == len(m) - 1:
            jm = sm[i]
            jx[i] = _sum(map(lambda a, b: a*b, x, m)) / sm[-1] * sp.sqrt(jm)
        else:
            jm = sm[i] * m[i + 1] / sm[i + 1]
            jx[i] = (_sum(map(lambda a, b: a*b, x[:i+1], m[:i+1]))
                / _sum(m[:i+1]) - x[i+1]) * sp.sqrt(jm)
    return jx

def toJacobiCoord(expr,x,m,jx):
    #there must be some nice way to make variables without names in sympy...
    uni='I_am_trying_to_break_things_if_I_choose_this_name_outside_of_toJacobiCoord'
    #solve one-D problem, use same solution many times
    x1=sp.symbols(['x'+uni+str(i) for i in range(len(x))])
    jx1=sp.symbols(['jxasdf'+uni+str(i) for i in range(len(x))])
    rep1 = sp.solve(map(lambda a, b: a-b, jx1, jacobiCoord(x1,m)),x1)
    #make extended dictonary for mD
    rep=dict()
    for i in rep1:
        for di in range(len(x[0])):
            toMultD=dict()
            for ii in range(len(x)):
                toMultD[x1[ii]]=x[ii][di]
                toMultD[jx1[ii]]=jx[ii][di]
            rep[i.subs(toMultD)]=rep1[i].subs(toMultD)
    return [e.subs(rep) for e in expr]
'''
n=5
d=3
#m=sp.symbols(['m'+str(i) for i in range(n)])
m=[sp.symbols('m')]*n
x=[[sp.symbols(['x'+str(i)+','+str(j)]) for j in range(d)] for i in range(n)]

jx=[[sp.symbols(['jx'+str(i)+','+str(j)]) for j in range(d)] for i in range(n)]
for i in range(n):
    for j in range(i):
        print(str(i)+", "+str(j)+":")
        print(sum( (x[i][di]-x[j][di])**2 for di in range(d) )  )
        print(toJacobiCoord([sum( (x[i][di]-x[j][di])**2 for di in range(d) )]
        ,x,m,jx))'''

'''print(sum(map(lambda a: a * a, jacobiCoord(x, m))))
print("asdf")
print(sp.simplify(sum(map(lambda a: a * a, jacobiCoord(x,
    m)))))'''

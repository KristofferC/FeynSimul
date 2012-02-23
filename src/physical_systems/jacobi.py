import sympy as sp

def jacobiCoord(x,m):
    jx=range(len(m))
    _sum = lambda a: reduce(lambda b,c:b+c,a)
    sm=[sum(m[:i+1]) for i in range(len(m))]
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

def toJacobiCoord(expr,x,m):
    jx = sp.symbols(['jx' + str(i) for i in range(len(m))])
    rep = sp.solve(map(lambda a, b: a-b, jx, jacobiCoord(x,m)),x)
    return [e.subs(rep) for e in expr]
'''
n=15
#m=sp.symbols(['m'+str(i) for i in range(n)])
m=[sp.symbols('m')]*n
x=sp.symbols(['x'+str(i) for i in range(n)])

jx = sp.symbols(['jx' + str(i) for i in range(len(m))])
rep = sp.solve(map(lambda a, b: a-b, jx, jacobiCoord(x,m)),x)
for i in range(n):
    for j in range(i):
        print(str(i)+", "+str(j)+":")
        print(((x[i]-x[j])**2).subs(rep).simplify())'''

'''print(sum(map(lambda a: a * a, jacobiCoord(x, m))))
print("asdf")
print(sp.simplify(sum(map(lambda a: a * a, jacobiCoord(x,
    m)))))'''

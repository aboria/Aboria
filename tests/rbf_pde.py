from __future__ import division
from sympy import MatrixSymbol, Matrix,diff,init_printing,symbols,sqrt,pprint,simplify


x = MatrixSymbol('x', 1, 2)

x = symbols('x')
y = symbols('y')
c = symbols('c')
g = sqrt(x**2 + y**2 + c**2)
gx = diff(g,x)
gy = diff(g,y)
gxx = diff(gx,x)
gyy = diff(gy,y)

print '---------------------------------------------'
pprint(g)
print '---------------------------------------------'
pprint(gx)
print '---------------------------------------------'
pprint(gy)
print '---------------------------------------------'
pprint(gxx)
print '---------------------------------------------'
pprint(gyy)
print '---------------------------------------------'
pprint(simplify(gxx+gyy))

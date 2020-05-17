from sympy import *

x = Symbol("x", real=True)
ω = Symbol("omega", real=True)

f = exp( -10*sin(ω*x)**2 )
print( diff(f,x,2) )
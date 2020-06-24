from sympy import *

A = symbols("A")
x_0 = symbols("x_0")
b = symbols("b")
c = symbols("c")
x = symbols("x")
Q = symbols("Q")
ρ = symbols("rho")
r_s = symbols("r_s")

def X(x):
    return x**2 + b*x + c

ε_c = A/2*( log(x/X(x)) + 2*b/Q*atan(Q/(2*x + b)) - b*x_0/X(x_0)*(
    log( (x-x_0)**2/X(x_0) ) + 2*(b + 2*x_0)/Q*atan(Q/(2*x + b))
    )
)

pprint(ε_c)

# Parameters
dict_subs_1 = {
    x: sqrt(r_s)
}

ε_c_subs1 = simplify(ε_c.subs(dict_subs_1))
pprint(ε_c_subs1)

pprint(diff(ε_c_subs1, r_s))

dict_subs_2 = {
    r_s: (4*pi*ρ/3)**Rational(1,3)
}

ε_c_subs2 = simplify(ε_c_subs1.subs(dict_subs_2))
pprint(ε_c_subs2)

pprint(diff(ε_c_subs2,ρ))

dict_subs_3 = {
    A: 0.0310907,
    x_0: -0.10498,
    b: 3.72744,
    c: 12.9352
}
ε_c_subs3 = simplify(ε_c_subs2.subs(dict_subs_3))
#pprint(ε_c_subs3)
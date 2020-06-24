from sympy import *

A = symbols("A", real=True)
x_0 = symbols("x_0", real=True)
b = symbols("b", real=True)
c = symbols("c", real=True)
x = symbols("x", real=True, positive=True)
Q = symbols("Q", real=True)
ρ = symbols("rho", real=True)
r_s = symbols("r_s", real=True, positive=True)

def X(x):
    return x**2 + b*x + c

ε_c = A*( log(x*x/X(x)) + 2*b/Q*atan(Q/(2*x + b)) - b*x_0/X(x_0)*(
    log( (x-x_0)**2/X(x) ) + 2*(b + 2*x_0)/Q*atan(Q/(2*x + b))
    )
)

#ε_c = A*( log(x*x/X(x)) + 2*b/Q*atan(Q/(2*x + b)) - (b*x_0)/X(x_0)*(log((x-x_0)*(x-x_0)/X(x)) + 2*(2*x_0+b)/Q*#atan(Q/(2*x + b)) ))

dict_subs_2 = {
    Q: sqrt(4*c - b**2)
}


dict_subs_3 = {
    A: 0.0310907,
    x_0: -0.10498,
    b: 3.72744,
    c: 12.9352,
}

pprint(ε_c.subs(dict_subs_3))

dict_subs_4 = {
    x: sqrt(r_s)
}

expr1 = ε_c.subs(dict_subs_2)
pprint("\n expr1 in x subs Q \n")
pprint(expr1)

expr1 = ε_c.subs(dict_subs_2).subs(dict_subs_3).subs(dict_subs_4) 
pprint("\n expr1 in r_s \n")
pprint(expr1)


"""
import numpy as np
import matplotlib.pyplot as plt
plt.clf()
Npts_plot = 400
r_s_num = np.linspace(0.001,100,Npts_plot)
eps_c_num = np.zeros(Npts_plot)
for ip in range(Npts_plot):
    eps_c_num[ip] = expr1.subs({r_s: r_s_num[ip]})
plt.plot(r_s_num, eps_c_num)
plt.savefig("IMG_eps_c_VWN.pdf")
"""

#d_V_c = diff(ε_c.subs(dict_subs_3), x)
d_V_c = diff(ε_c.subs(dict_subs_2), x)
print("\nd_V_c = \n")
pprint(simplify(d_V_c))

print("\nd_V_c subs Q = \n")
pprint( simplify(d_V_c.subs( {Q: sqrt(4*c - b**2)} ) ) )

expr1 = A*( 2/x - (2*x+b)/X(x) - 4*b/(Q*Q + (2*x+b)*(2*x+b)) - (b*x_0)/X(x_0)*(
    2/(x-x_0) - (2*x+b)/X(x) - 4*(2*x_0+b) / (Q*Q + (2*x+b) * (2*x+b)) ) )
#expr1 = expr1.subs(dict_subs_3)
print("\nexpr1 = \n")
pprint(simplify(expr1))

print("\nexpr1 subs Q = \n")
pprint(simplify(expr1.subs( {Q: sqrt(4*c - b**2)} ) ) )

"""
plt.clf()
Npts_plot = 400
r_s_num = np.linspace(0.001,1,Npts_plot)
d_V_c_num = np.zeros(Npts_plot)
d_V_c_num2 = np.zeros(Npts_plot)
for ip in range(Npts_plot):
    d_V_c_num[ip] = expr1.subs({r_s: r_s_num[ip]})
    d_V_c_num2[ip] = d_V_c.subs({r_s: r_s_num[ip]})
plt.plot(r_s_num, d_V_c_num, label="v1")
plt.plot(r_s_num, d_V_c_num2, label="v2")
plt.legend()
plt.grid()
plt.savefig("IMG_d_V_c_VWN.pdf")

diff_V = expr1 - d_V_c
for r in [0.01, 0.02, 0.4, 0.6, 1.0]:
    print(diff_V.subs({r_s: r}))

print("expr1 = ")
pprint(simplify(expr1))
print("d_V_c = ")
pprint(simplify(d_V_c))

#d_V_c = simplify(d_V_c)
#pprint(simplify(ε_c*d_V_c))
#for term1 in d_V_c.args:
#    pprint(term1)

#pprint( simplify( d_V_c.args[-1].expand() ) )

"""

"""
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
"""
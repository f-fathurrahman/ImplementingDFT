from sympy import *
#from integer_numeric import *
from integer_symbols import *

X_FACTOR_C = THREE/EIGHT*(THREE/pi)**(ONE/THREE) * 4**(TWO/THREE)

# Dimension = 3
RS_FACTOR = (THREE/(FOUR*pi))**(ONE/THREE)
LDA_X_FACTOR = -X_FACTOR_C

params_a_alpha = 1
lda_x_ax = -params_a_alpha*RS_FACTOR*X_FACTOR_C/2**(FOUR/THREE)

def f_lda_x(r_s, z):
    return lda_x_ax*( (1 + z)**(FOUR/THREE) + (1 - z)**(FOUR/THREE) )/r_s

def f(r_s, z):
    return f_lda_x(r_s, z)


ζ = symbols("zeta")
ρ = symbols("rho")

r_s = (THREE/FOUR/pi)**(ONE/THREE)*ρ**(-ONE/THREE)
#r_s = symbols("r_s")

print()
print("Spin pol")
pprint(f(r_s, ζ))

print()
print("Non-spin pol")
ζ = 0
pprint(f(r_s, ζ))

E_x = f(r_s, ζ)

print("Derivative w.r.t ρ:")
pprint(diff(E_x, ρ))

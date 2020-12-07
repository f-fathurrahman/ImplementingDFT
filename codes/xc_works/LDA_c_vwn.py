from sympy import *
from integer_numeric import *
#from integer_symbols import *

# from util
# See Eq. (9) of Perdew1992_13244
def f_zeta(z):
    return ((1 + z)**(4/3) + (1 - z)**(4/3) - 2)/(2**(4/3) - 2)

# zero elements are set to zero (not used)
A_vwn  = [0.0,  0.0310907, 0.01554535, -1/(6*pi**2)]
b_vwn  = [0.0,  3.72744,   7.06042,    1.13107  ]
c_vwn  = [0.0, 12.9352,   18.0578,    13.0045   ]
x0_vwn = [0.0, -0.10498,  -0.32500,   -0.0047584]

def Q_vwn(b, c):
    return sqrt(4*c - b**2)

def f1_vwn(b, c):
    return 2*b/Q_vwn(b, c)

def f2_vwn(b, c, x0):
    return b*x0/(x0**2 + b*x0 + c)

def f3_vwn(b, c, x0):
    return 2*(2*x0 + b)/Q_vwn(b, c)

fpp_vwn = 4/(9*(2**(1/3) - 1))

def fx_vwn(b, c, rs):
    return rs + b*sqrt(rs) + c

def DMC(rs, z):
    return f_aux(A_vwn[2], b_vwn[2], c_vwn[2], x0_vwn[2], rs) - \
        f_aux(A_vwn[1], b_vwn[1], c_vwn[1], x0_vwn[1], rs)

def f_aux(A, b, c, x0, rs):
    return A*( log(rs/fx_vwn(b, c, rs)) +
        (f1_vwn(b, c) - f2_vwn(b, c, x0)*f3_vwn(b, c, x0))*atan(Q_vwn(b, c)/(2*sqrt(rs) + b)) -
         f2_vwn(b, c, x0)*log((sqrt(rs) - x0)**2/fx_vwn(b, c, rs))
    )

def f_vwn(rs, z):
    return f_aux(A_vwn[1], b_vwn[1], c_vwn[1], x0_vwn[1], rs) + \
        f_aux(A_vwn[3], b_vwn[3], c_vwn[3], x0_vwn[3], rs)*f_zeta(z)*(1 - z**4)/fpp_vwn + \
        DMC(rs, z)*f_zeta(z)*z**4

def f(rs, z):
    return f_vwn(rs, z)


ζ = symbols("zeta", real=True)
ρ = symbols("rho", real=True)

#r_s = (THREE/FOUR/pi)**(ONE/THREE)*ρ**(-ONE/THREE)
r_s = symbols("r_s", real=True)

print()
print("Spin pol")
ec_spin = simplify(f(r_s, ζ))
pprint(ec_spin)

print()
print("Non-spin pol")
ζ = 0
ec_nonspin = simplify(f(r_s, ζ))
pprint(ec_nonspin)

from sympy.utilities.codegen import codegen

code1 = codegen( ("ec_spin", ec_spin), language="julia")
print(code1[0][1])

code1 = codegen( ("ec_nonspin", ec_nonspin), language="julia")
print(code1[0][1])


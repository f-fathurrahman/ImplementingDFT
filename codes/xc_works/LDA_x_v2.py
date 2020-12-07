from sympy import *
#from integer_numeric import *
from integer_symbols import *

Rhoe = symbols("ρ")

third = ONE/THREE
pi34 = (THREE/(FOUR*pi))**(ONE/THREE)
f = -NINE/EIGHT*(THREE/TWO/pi)**(TWO/THREE)
alpha = TWO/THREE

r_s = pi34/Rhoe**third
#r_s = symbols("r_s")

ex = f * alpha / r_s

pprint(ex)
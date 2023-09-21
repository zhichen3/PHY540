'''
Solve for temperature:
3678 = \sum E_i exp(-E_i/(kB*T)) / \sum exp(-E_i/(kB*T))
'''

from scipy.optimize import fsolve
import numpy as np
from sympy import *
init_printing()

# Energies
E_i = np.array([251.31, 299.42, 417.20, 829.49, 911.66,
                1042.48, 1124.40, 1194.48, 1291.82, 1310.12,
                1427.20, 1481.76, 1513.84, 1531.49, 1560.50,
                2982.59, 3006.07, 3058.31, 3129.27, 3134.26,
                3752.52])

# Create symbols for sympy
b = symbols('b')
E = symbols('E')
lnQ = symbols('lnQ')
N = symbols('N')
i = symbols('i')

# Sum(ln(Q_i)) which is the same as ln(Prod(Q_i))
lnQ = Sum(ln(1.0/(1 - exp(-b*Indexed('E', i)))), (i, 0, len(E_i)-1))

# Differentiate w.r.t beta
dlnQdBeta = diff(lnQ, b)

# Equation we want to solve for is:
# dlnQdBeta + 3678 = 0
f = dlnQdBeta + 3678.0

# Print out the equation we want to solve for.
pprint(f)

# Lambdify sympy expressions so we can use numpy.fsolve to solve for beta
func = lambdify((b, E), f)

# Solve for beta
sol = fsolve(func, 0.001, args=(E_i))

# Convert beta solution to temperature
T = 1.0 / (sol * 0.695)

# Print out solution
print(f"beta is {sol[0]}")
print(f"T is {T[0]}")


import sympy as sp

# Define the variables
x, y, x0, y0, sigma = sp.symbols('x y x0 y0 sigma')

# Define the functions
Lambda_x = sp.erf((x + 1/2 - x0)/(sp.sqrt(2)*sigma)) - sp.erf((x - 1/2 - x0)/(sp.sqrt(2)*sigma))
Lambda_y = sp.erf((y + 1/2 - y0)/(sp.sqrt(2)*sigma)) - sp.erf((y - 1/2 - y0)/(sp.sqrt(2)*sigma))
L = Lambda_x * Lambda_y

# Compute the common factors
j_x0 = L.diff(x0)
j_y0 = L.diff(y0)
j_sigma = L.diff(sigma)

# Compute the Hessian using the common factors
J = sp.Matrix([j_x0,j_y0,j_sigma])

# Convert to numerical function
J_func = sp.lambdify((x, y, x0, y0, sigma), J, "numpy")

# Generate code for the Jacobian function
code = """import numpy as np\n
from numpy import sqrt\n
from numpy import exp\n
from numpy import pi\n
from scipy.special import erf\n\n
"""

code += "def jacobian1(x, y, x0, y0, sigma):\n"
code += f"    j_x0 = {j_x0}\n"
code += f"    j_y0 = {j_y0}\n"
code += f"    j_sigma = {j_sigma}\n"
code += "    jac = np.array([j_x0, j_y0, j_sigma], dtype=np.float64)\n"
code += "    return jac\n"

# Save code to file
with open("jac1.py", "w") as f:
    f.write(code)


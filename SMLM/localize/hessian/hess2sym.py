import sympy as sp

# Define the variables
x, y, x0, y0, sigma = sp.symbols('x y x0 y0 sigma')

# Define the functions
Lambda_x = sp.erf((x + 1/2 - x0)/(sp.sqrt(2)*sigma)) - sp.erf((x - 1/2 - x0)/(sp.sqrt(2)*sigma))
Lambda_y = sp.erf((y + 1/2 - y0)/(sp.sqrt(2)*sigma)) - sp.erf((y - 1/2 - y0)/(sp.sqrt(2)*sigma))
L = Lambda_x * Lambda_y

# Compute the common factors
h_xx = L.diff(x0).diff(x0)
h_xy = L.diff(x0).diff(y0)
h_xs = L.diff(x0).diff(sigma)
h_yx = L.diff(y0).diff(x0)
h_yy = L.diff(y0).diff(y0)
h_ys = L.diff(y0).diff(sigma)
h_sx = L.diff(sigma).diff(x0)
h_sy = L.diff(sigma).diff(y0)
h_ss = L.diff(sigma).diff(sigma)

# Compute the Hessian using the common factors
H = sp.Matrix([[h_xx,h_xy,h_xs], [h_yx,h_yy,h_ys], [h_sx,h_sy,h_ss]])

# Convert to numerical function
H_func = sp.lambdify((x, y, x0, y0, sigma), H, "numpy")

# Generate code for the Hessian function
code = """import numpy as np\n
from numpy import sqrt\n
from numpy import exp\n
from numpy import pi\n
from scipy.special import erf\n\n
"""

code += "def hessian2(x, y, x0, y0, sigma):\n"
code += f"    h_xx = {h_xx}\n"
code += f"    h_xy = {h_xy}\n"
code += f"    h_xs = {h_xs}\n"
code += f"    h_yx = {h_yx}\n"
code += f"    h_yy = {h_yy}\n"
code += f"    h_ys = {h_ys}\n"
code += f"    h_sx = {h_sx}\n"
code += f"    h_sy = {h_sy}\n"
code += f"    h_ss = {h_ss}\n"
code += "    hessian = np.array([[h_xx,h_xy,h_xs], [h_yx,h_yy,h_ys], [h_sx,h_sy,h_ss]], dtype=np.float64)\n"
code += "    return hessian\n"

# Save code to file
with open("hess2.py", "w") as f:
    f.write(code)


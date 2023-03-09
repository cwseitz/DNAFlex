import sympy as sp

# Define the variables
x, y, x0, y0, sigma, N0 = sp.symbols('x y x0 y0 sigma N0')

# Define the functions
Lambda_x = sp.erf((x + 1/2 - x0)/(sp.sqrt(2)*sigma)) - sp.erf((x - 1/2 - x0)/(sp.sqrt(2)*sigma))
Lambda_y = sp.erf((y + 1/2 - y0)/(sp.sqrt(2)*sigma)) - sp.erf((y - 1/2 - y0)/(sp.sqrt(2)*sigma))
L = 0.25*Lambda_x*Lambda_y

# Compute elements of the Hessian
h_xx = L.diff(x0).diff(x0)
h_xy = L.diff(x0).diff(y0)
h_xs = L.diff(x0).diff(sigma)
h_xN = L.diff(x0).diff(N0)

h_yx = L.diff(y0).diff(x0)
h_yy = L.diff(y0).diff(y0)
h_ys = L.diff(y0).diff(sigma)
h_yN = L.diff(y0).diff(N0)

h_sx = L.diff(sigma).diff(x0)
h_sy = L.diff(sigma).diff(y0)
h_ss = L.diff(sigma).diff(sigma)
h_sN = L.diff(sigma).diff(N0)

h_Nx = L.diff(N0).diff(x0)
h_Ny = L.diff(N0).diff(y0)
h_Ns = L.diff(N0).diff(sigma)
h_NN = L.diff(N0).diff(N0)

H = sp.Matrix([[h_xx,h_xy,h_xs,h_xN], [h_yx,h_yy,h_ys,h_yN], [h_sx,h_sy,h_ss,h_sN],[h_Nx,h_Ny,h_Ns,h_NN]])

# Convert to numerical function
H_func = sp.lambdify((x, y, x0, y0, sigma, N0), H, "numpy")

# Generate code for the Hessian function
code = """import numpy as np\n
from numpy import sqrt\n
from numpy import exp\n
from numpy import pi\n
from scipy.special import erf\n\n
"""

code += "def hessian2(x, y, x0, y0, sigma, N0):\n"
code += f"    h_xx = {h_xx}\n"
code += f"    h_xy = {h_xy}\n"
code += f"    h_xs = {h_xs}\n"
code += f"    h_xN = {h_xN}\n"
code += f"    h_yx = {h_yx}\n"
code += f"    h_yy = {h_yy}\n"
code += f"    h_ys = {h_ys}\n"
code += f"    h_yN = {h_yN}\n"
code += f"    h_sx = {h_sx}\n"
code += f"    h_sy = {h_sy}\n"
code += f"    h_ss = {h_ss}\n"
code += f"    h_sN = {h_sN}\n"
code += f"    h_Nx = {h_Nx}\n"
code += f"    h_Ny = {h_Ny}\n"
code += f"    h_Ns = {h_Ns}\n"
code += f"    h_NN = {h_NN}\n"
code += "    H = np.array([[h_xx,h_xy,h_xs,h_xN], [h_yx,h_yy,h_ys,h_yN], [h_sx,h_sy,h_ss,h_sN],[h_Nx,h_Ny,h_Ns,h_NN]], dtype=np.float64)\n"
code += "    return hessian\n"

# Save code to file
with open("hess2.py", "w") as f:
    f.write(code)


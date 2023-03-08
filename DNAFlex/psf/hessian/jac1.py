import numpy as np
from numpy import sqrt
from numpy import exp
from numpy import pi
from scipy.special import erf

def jacobian1(x, y, x0, y0, sigma):
    j_x0 = (-sqrt(2)*exp(-(x - x0 + 0.5)**2/(2*sigma**2))/(sqrt(pi)*sigma) + sqrt(2)*exp(-(x - x0 - 0.5)**2/(2*sigma**2))/(sqrt(pi)*sigma))*(-erf(sqrt(2)*(y - y0 - 0.5)/(2*sigma)) + erf(sqrt(2)*(y - y0 + 0.5)/(2*sigma)))
    j_y0 = (-sqrt(2)*exp(-(y - y0 + 0.5)**2/(2*sigma**2))/(sqrt(pi)*sigma) + sqrt(2)*exp(-(y - y0 - 0.5)**2/(2*sigma**2))/(sqrt(pi)*sigma))*(-erf(sqrt(2)*(x - x0 - 0.5)/(2*sigma)) + erf(sqrt(2)*(x - x0 + 0.5)/(2*sigma)))
    j_sigma = (sqrt(2)*(x - x0 - 0.5)*exp(-(x - x0 - 0.5)**2/(2*sigma**2))/(sqrt(pi)*sigma**2) - sqrt(2)*(x - x0 + 0.5)*exp(-(x - x0 + 0.5)**2/(2*sigma**2))/(sqrt(pi)*sigma**2))*(-erf(sqrt(2)*(y - y0 - 0.5)/(2*sigma)) + erf(sqrt(2)*(y - y0 + 0.5)/(2*sigma))) + (sqrt(2)*(y - y0 - 0.5)*exp(-(y - y0 - 0.5)**2/(2*sigma**2))/(sqrt(pi)*sigma**2) - sqrt(2)*(y - y0 + 0.5)*exp(-(y - y0 + 0.5)**2/(2*sigma**2))/(sqrt(pi)*sigma**2))*(-erf(sqrt(2)*(x - x0 - 0.5)/(2*sigma)) + erf(sqrt(2)*(x - x0 + 0.5)/(2*sigma)))
    jac = np.array([j_x0, j_y0, j_sigma], dtype=np.float64)
    return jac

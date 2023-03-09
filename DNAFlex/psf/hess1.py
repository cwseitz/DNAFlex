import numpy as np
from scipy.special import erf


def hessian1(counts,x,y,x0,y0,sigma):
    alpha = np.sqrt(2)*sigma
    lamdx = 0.5*erf((x+0.5-x0)/alpha) - erf((x-0.5-x0)/alpha)
    lamdy = 0.5*erf((y+0.5-y0)/alpha) - erf((y-0.5-y0)/alpha)
    lamd = lamdx*lamdy
    hess = np.diag(-1*counts/(lamd**2))
    return hess

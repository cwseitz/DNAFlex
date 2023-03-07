import numpy as np
from scipy.special import erf

def jacobian2(counts,x,y,x0,y0,sigma):
    alpha = np.sqrt(2)*sigma
    lamdx = erf((x+0.5-x0)/alpha) - erf((x-0.5-x0)/alpha)
    lamdy = erf((y+0.5-y0)/alpha) - erf((y-0.5-y0)/alpha)
    lamd = lamdx*lamdy
    jac = counts/lamd - 1
    return jac

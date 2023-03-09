import numpy as np
from scipy.special import erf

def jacobian2(counts,x,y,x0,y0,sigma,N0,texp,eta):
    alpha = np.sqrt(2)*sigma
    lamdx = 0.5*erf((x+0.5-x0)/alpha) - erf((x-0.5-x0)/alpha)
    lamdy = 0.5*erf((y+0.5-y0)/alpha) - erf((y-0.5-y0)/alpha)
    lamd = lamdx*lamdy
    I0 = eta*texp*N0
    jac = counts/lamd - I0
    return jac

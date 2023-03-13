import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

from .hess1 import hessian1
from .hess2 import hessian2 
from .jac1 import jacobian1
from .jac2 import jacobian2

def hessian_analytical(theta,counts,eta,texp):
    lx, ly = counts.shape
    x0,y0,sigma,N0 = theta
    X,Y = np.meshgrid(np.arange(0,lx),np.arange(0,ly))
    X = X.ravel(); Y = Y.ravel(); counts = counts.ravel()
    J1 = jacobian1(X,Y,x0,y0,sigma,N0)
    H1 = hessian1(counts,X,Y,x0,y0,sigma,N0)
    J2 = jacobian2(counts,X,Y,x0,y0,sigma,N0,eta,texp)
    H2 = hessian2(X,Y,x0,y0,sigma,N0)
    A = J1 @ H1 @ J1.T 
    B = np.sum(H2*J2[np.newaxis, np.newaxis, :],axis=-1)
    return A + B

#def negloglike(theta,counts,eta,texp):
#    lx, ly = counts.shape
#    x0,y0,sigma,N0 = theta
#    alpha = np.sqrt(2)*sigma
#    X,Y = np.meshgrid(np.arange(0,lx),np.arange(0,ly))
#    X = X.ravel(); Y = Y.ravel(); counts = counts.ravel()
#    lamdx = 0.5*(erf((X+0.5-x0)/alpha) - erf((X-0.5-x0)/alpha))
#    lamdy = 0.5*(erf((Y+0.5-y0)/alpha) - erf((Y-0.5-y0)/alpha))
#    I0 = eta*N0*texp
#    mu = I0*lamdx*lamdy + 1e-8
#    counts = counts + 1e-8
#    stirling = counts*np.log(counts) - counts
#    ll = counts*np.log(mu) - stirling - mu
#    ll = np.sum(ll)
#    return -1*ll

def negloglike(theta,counts,eta,texp):
    lx, ly = counts.shape
    x0,y0,sigma,N0,B0 = theta
    alpha = np.sqrt(2)*sigma
    X,Y = np.meshgrid(np.arange(0,lx),np.arange(0,ly))
    X = X.ravel(); Y = Y.ravel(); counts = counts.ravel()
    lamdx = 0.5*(erf((X+0.5-x0)/alpha) - erf((X-0.5-x0)/alpha))
    lamdy = 0.5*(erf((Y+0.5-y0)/alpha) - erf((Y-0.5-y0)/alpha))
    I0 = eta*N0*texp
    B = eta*B0*texp
    mu = I0*lamdx*lamdy + B + 1e-8
    counts = counts + 1e-8
    stirling = counts*np.log(counts) - counts
    ll = counts*np.log(mu) - stirling - mu
    ll = np.sum(ll)
    return -1*ll




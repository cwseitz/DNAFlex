import autograd.numpy as np
import matplotlib.pyplot as plt
from autograd import grad, jacobian, hessian
from autograd.scipy.stats import norm, multivariate_normal
from autograd.scipy.special import erf
from scipy.optimize import minimize
from scipy.special import factorial
np.set_printoptions(suppress=True)
np.random.seed(10)

from .hess1 import hessian1
from .hess2 import hessian2 
from .jac1 import jacobian1
from .jac2 import jacobian2
import matplotlib.pyplot as plt

def hessian_analytical(theta, counts):
    lx, ly = counts.shape
    X,Y = np.meshgrid(np.arange(0,lx),np.arange(0,ly))
    X = X.ravel(); Y = Y.ravel(); counts = counts.ravel()
    J1 = jacobian1(X,Y,*theta)
    H1 = hessian1(counts,X,Y,*theta)
    J2 = jacobian2(counts,X,Y,*theta)
    H2 = hessian2(X,Y,*theta)
    A = J1 @ H1 @ J1.T 
    B = np.sum(H2*J2[np.newaxis, np.newaxis, :],axis=-1)
    return A + B

#def hessian_autograd(theta):
#    def loglike(theta):
#        x0,y0,sigma = theta
#        alpha = np.sqrt(2)*sigma
#        lamdx = erf((X+0.5-x0)/alpha) - erf((X-0.5-x0)/alpha)
#        lamdy = erf((Y+0.5-y0)/alpha) - erf((Y-0.5-y0)/alpha)
#        lamd = I0*lamdx*lamdy
#        ll = counts*np.log(lamd) - np.log(factorial(counts)) - lamd
#        ll = np.sum(ll)
#        return ll
#    hessian_ = hessian(loglike)
#    hess = hessian_(theta)
#    return hess




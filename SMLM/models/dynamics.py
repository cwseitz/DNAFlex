import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from copy import deepcopy
from matplotlib import cm
from scipy.stats import multivariate_normal

class BrownianMotion:
    def __init__(self, N, T, dt, tau, cov, dx=0.1, x_max=1, x0=None, trials=1000, dtype=np.float32):

        """
        Integrate N-dimensional Langevin dynamics for stationary delta-correlated Gaussian noise with zero mean (Brownian motion)
        """
        
        #Params
        self.N = N
        self.nsteps = int(round(T/dt))
        self.dt = dt
        self.dx = dx
        self.x0 = x0
        self.x_max = x_max
        self.nx = int(round(2*self.x_max/self.dx))
        self.trials = trials
        self.mu =  np.zeros((N,))
        self.cov = cov

        self.X = np.zeros((self.trials, self.nsteps, self.N))
        if self.x0 != None:
             self.X[:,:,0] = self.x0

    def simulate(self):
        f = multivariate_normal(mean=self.mu,cov=self.cov)
        dW = f.rvs(size=(self.trials,self.nsteps))
        for i in range(self.trials):
            for j in range(1,self.nsteps):
                self.X[i,j,:] = self.X[i,j-1,:] + dW[i,j,:]*np.sqrt(self.dt)
                           
class Poisson:

    """

    Generate an ensemble of Poisson spike trains from a vector of rates

    Parameters
    ----------
    rates : ndarray
        A matrix - each value providing the firing rate for a single input unit
        By default, generate an ensemble of homogeneous poisson processes
        with rate 20Hz

    Returns
    -------
    spikes : 2d ndarray
        a binary matrix containing the spiking patterns for the ensemble
        (i, j) --> (unit, nsteps)

    """

    def __init__(self, T, dt, N, trials=1, rates=None, random_select=None):

        self.T = T
        self.dt = dt
        self.N = N
        self.nsteps = 1 + int(round(T/dt))
        self.random_select = random_select
        self.trials = trials

        if rates is None:
            self.r0 = 20 #default rate (Hz)
            rates = self.r0*np.ones((self.N, self.trials, self.nsteps))

        self.rates = rates
        self.run_generator()

    def run_generator(self):

        self.r = self.rates*self.dt
        self.x = np.random.uniform(0,1,size=(self.N,self.trials,self.nsteps))
        self.spikes = np.array(self.x < self.r, dtype=np.int32)

        if self.random_select != None:
            rng = default_rng()
            x = rng.choice(self.N, size=self.random_select, replace=False)
            self.spikes[x,:,:] = 0

    def to_currents(self, J):
        self.currents = np.einsum('ij,jhk->ihk', J, self.spikes)
        return self.currents

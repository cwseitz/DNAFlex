import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, r_
from scipy import optimize
from .model import *

class Gaussian2D(Model):

    def __init__(self):
        super().__init__()

    def gaussian(self, A, x0, y0, sigma):
        sigma = float(sigma)
        return lambda x,y: A*np.exp(
                    -(((x-x0)/sigma)**2+((y-y0)/sigma)**2)/2)

    def moments(self, data):
        total = data.sum()
        X, Y = np.indices(data.shape)
        x = (X*data).sum()/total
        y = (Y*data).sum()/total
        col = data[:, int(y)]
        sigma = np.sqrt(np.abs((np.arange(col.size)-x)**2*col).sum()/col.sum())
        row = data[int(x), :]
        A = data.max()
        return A, x, y, sigma

    def fit(self, data):
        theta = self.moments(data)
        errorfunction = lambda p: np.ravel(self.gaussian(*p)(*np.indices(data.shape)) -
                                     data)
        p, success = optimize.leastsq(errorfunction, theta)
        return p

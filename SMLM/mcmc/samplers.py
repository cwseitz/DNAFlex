import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal

class Sampler:
    def __init__(self):
        pass
        
class MetropolisHastings(Sampler):
    def __init__(self,mu,cov,likelihood,prior):
        super().__init__()
        self.mu = mu
        self.cov = cov
        self.likelihood = likelihood
        self.prior = prior
    def sample(self,nsamples,x0):
        samples = []
        x = x0
        for i in range(nsamples):
            dx = np.random.multivariate_normal(mean=self.mu, cov=self.cov)
            x_new = x + dx
            like1 = self.likelihood.eval(x_new)
            prior1 = self.prior.eval(x_new[0],x_new[1])
            like2 = self.likelihood.eval(x)
            prior2 = self.prior.eval(x[0],x[1])
            a1 = (like1*prior1)/(like2*prior2)
            a2 = 1; a = a1*a2
            if a >= 1:
                x = x_new
                u = None
            else:
                u = np.random.uniform(0,1)
                if u <= a:
                    x = x_new
            samples.append(x)
            print(f'Iteration: {i}, m={x[0]}, b={x[1]}')
        return np.array(samples)
        
           

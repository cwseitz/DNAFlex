from scipy.stats import poisson, norm
from smlm.models import PoissonNormal
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(0,30)
fig, ax = plt.subplots()
mu_psn = 10
mu_norm, sigma_norm = 10,0.1
pnorm = PoissonNormal(mu_norm,mu_psn,sigma_norm)
X = pnorm.sample(20000)
y = pnorm.get_pmf(x)
plt.hist(X,density=True,color='red')
plt.plot(x,y,color='blue')

#Approximate PoissonNormal as Normal
var = sigma_norm**2 + mu_psn
g = norm(loc=mu_norm+mu_psn,scale=np.sqrt(var))
gpdf = g.pdf(x)
plt.plot(x,gpdf,linestyle='--',color='cyan')
plt.show()

plt.plot(np.cumsum(y))
plt.plot(np.cumsum(gpdf))
plt.show()

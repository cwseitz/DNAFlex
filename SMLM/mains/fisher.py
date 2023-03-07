from smlm.models import NoisyGaussianPSFAnalytic, NoisyGaussianPSFMLE
import matplotlib.pyplot as plt
import numpy as np

L = 19
omat = np.ones((L,L))
gain0 = 2.2
rmu0 = 100
rvar0 = 700
gain = gain0*omat #ADU/e-
rmu = rmu0*omat #ADU
rvar = rvar0*omat #ADU^2
pixel_size = 108.3 #nm
sigma = 0.22*640/1.4 #zhang 2007
sigma = sigma = sigma/pixel_size
alpha = 1/(np.sqrt(2)*sigma)
lam0 = 10000 #cps
texp = 1 #seconds
eta = 1
b = 0
I0 = lam0*texp #expected number of photons
x0,y0 = (10,10)

npsf = NoisyGaussianPSFAnalytic(x0,y0,b,L,I0,alpha,eta,gain,rmu,rvar)
H = npsf.generate(plot=True)


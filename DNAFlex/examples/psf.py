from DNAFlex.psf import GaussianPSF

x0, y0 = 5, 5
L = 10
I0 = 1000
sigma = 1.2
eta = 1.0
gain = 1.0
rmu = 0
rvar = 2
psf = GaussianPSF(x0,y0,L,I0,sigma,eta,gain,rmu,rvar)
psf.generate(plot=True)

import autograd.numpy as np
import matplotlib.pyplot as plt
from autograd import grad, jacobian, hessian
from autograd.scipy.stats import norm, multivariate_normal
from autograd.scipy.special import erf
from scipy.optimize import minimize
from scipy.special import factorial

class NoisyGaussianPSF:
    def __init__(self,x0,y0,b,L,I0,alpha,eta,gain,rmu,rvar,depth=16):
        frame = np.zeros((L,L))
        x = np.arange(0,L); y = np.arange(0,L)
        X,Y = np.meshgrid(x,y)
        intx = 0.5*(erf((X+0.5-x0)*alpha)-erf((X-0.5-x0)*alpha))
        inty = 0.5*(erf((Y+0.5-y0)*alpha)-erf((Y-0.5-y0)*alpha))
        self.mu = b + eta*I0*intx*inty
        self.pe = np.random.poisson(lam=self.mu) #photoelectrons
        self.X = gain*self.pe
        self.read_noise = np.random.normal(rmu,np.sqrt(rvar))
        self.X += self.read_noise
        max_adu = 2**depth - 1
        self.X = self.X.astype(np.int16)
        self.X[self.X > max_adu] = max_adu # models pixel saturation

    def plot_generated(self):
        fig, ax = plt.subplots(2,2,figsize=(4,3))
        im1 = ax[0,0].imshow(self.mu,cmap='gray')
        ax[0,0].set_xticks([]);ax[0,0].set_yticks([])
        plt.colorbar(im1, ax=ax[0,0], label=r'$\mu_{p}$')
        im2 = ax[0,1].imshow(self.pe,cmap='gray')
        ax[0,1].set_xticks([]);ax[0,1].set_yticks([])
        plt.colorbar(im2, ax=ax[0,1], label=r'$e^{-}$')
        im3 = ax[1,0].imshow(self.read_noise,cmap='gray')
        ax[1,0].set_xticks([]);ax[1,0].set_yticks([])
        plt.colorbar(im3, ax=ax[1,0], label=r'$\xi$ (ADU)')
        im4 = ax[1,1].imshow(self.X,cmap='gray')
        ax[1,1].set_xticks([]);ax[1,1].set_yticks([])
        plt.colorbar(im4, ax=ax[1,1], label=r'ADU')
        plt.tight_layout()
        plt.show()

class NoisyGaussianPSFAnalytic:
    def __init__(self,x0,y0,b,L,I0,alpha,eta,gain,rmu,rvar,depth=16):
        self.x0 = x0 #true x-coordinate
        self.y0 = y0 #true y-coordinate
        self.b = b #number of background photons
        self.L = L #number of pixels along an axis
        self.I0 = I0 #number of photons emitted
        self.alpha = alpha
        self.eta = eta #quantum efficiency
        self.gain = gain #ADU/e-
        self.rmu = rmu #readout noise mean
        self.rvar = rvar #readout noise variance
    def generate(self,depth=16,plot=False):
        frame = np.zeros((self.L,self.L))
        x = np.arange(0,self.L); y = np.arange(0,self.L)
        X,Y = np.meshgrid(x,y)
        mux = 0.5*(erf((X+0.5-self.x0)*self.alpha)-erf((X-0.5-self.x0)*self.alpha))
        muy = 0.5*(erf((Y+0.5-self.y0)*self.alpha)-erf((Y-0.5-self.y0)*self.alpha))
        mu = self.b + self.eta*self.I0*mux*muy
        pe = np.random.poisson(lam=mu) #photoelectrons
        X = self.gain*pe
        read_noise = np.random.normal(self.rmu,np.sqrt(self.rvar))
        X += read_noise
        max_adu = 2**depth - 1
        X = X.astype(np.int16)
        X[X > max_adu] = max_adu # models pixel saturation
        if plot:
            self.plot_generated(mu,pe,read_noise,X)
        return pe
    def plot_generated(self,mu,pe,read_noise,X):
        fig, ax = plt.subplots(2,2,figsize=(4,3))
        im1 = ax[0,0].imshow(mu,cmap='gray')
        ax[0,0].set_xticks([]);ax[0,0].set_yticks([])
        plt.colorbar(im1, ax=ax[0,0], label=r'$\mu_{p}$')
        im2 = ax[0,1].imshow(pe,cmap='gray')
        ax[0,1].set_xticks([]);ax[0,1].set_yticks([])
        plt.colorbar(im2, ax=ax[0,1], label=r'$e^{-}$')
        im3 = ax[1,0].imshow(read_noise,cmap='gray')
        ax[1,0].set_xticks([]);ax[1,0].set_yticks([])
        plt.colorbar(im3, ax=ax[1,0], label=r'$\xi$ (ADU)')
        im4 = ax[1,1].imshow(X,cmap='gray')
        ax[1,1].set_xticks([]);ax[1,1].set_yticks([])
        plt.colorbar(im4, ax=ax[1,1], label=r'ADU')
        plt.tight_layout()
        plt.show()
    def neg_loglike(self,H,theta):
        x0,y0 = theta
        x = np.arange(0,self.L)
        y = np.arange(0,self.L)
        X,Y = np.meshgrid(x,y)
        mux = 0.5*(erf((X+0.5-x0)*self.alpha)-erf((X-0.5-x0)*self.alpha))
        muy = 0.5*(erf((Y+0.5-y0)*self.alpha)-erf((Y-0.5-y0)*self.alpha))
        mu = self.I0*mux*muy #assume we know I0 for now
        nk = H.flatten()
        muk = mu.flatten()
        s = nk*np.log(muk)
        nll = np.sum(s)
        return nll
    def hessian(self,H,x0,y0):
        sigma = 1/(np.sqrt(2)*self.alpha) #just for simplicity
        x = np.arange(0,self.L); y = np.arange(0,self.L)
        X,Y = np.meshgrid(x,y)
        mux = 0.5*(erf((X+0.5-x0)*self.alpha)-erf((X-0.5-x0)*self.alpha))
        muy = 0.5*(erf((Y+0.5-y0)*self.alpha)-erf((Y-0.5-y0)*self.alpha))
        mu = self.I0*mux*muy #assume we know I0 for now
        Jm = np.zeros((self.L**2,2))
        Jm[:,0] = (muy*np.sqrt(2/(np.pi*sigma**2))*(self.gauss(X,x0+0.5,sigma)-self.gauss(X,x0-0.5,sigma))).flatten()
        Jm[:,1] = (mux*np.sqrt(2/(np.pi*sigma**2))*(self.gauss(Y,y0+0.5,sigma)-self.gauss(Y,y0-0.5,sigma))).flatten()
        Hm = np.zeros((self.L**2,2,2))
        Hm[:,0,0] = (muy*np.sqrt(2/np.pi)*((X-x0-0.5)*self.gauss(X,x0+0.5,sigma)-(X-x0+0.5)*self.gauss(X,x0-0.5,sigma))/(sigma**3)).flatten()
        Hm[:,0,1] = (Jm[:,0]*Jm[:,1]).flatten()
        Hm[:,1,0] = Hm[:,0,1]
        Hm[:,1,1] = (mux*np.sqrt(2/np.pi)*((Y-y0-0.5)*self.gauss(Y,y0+0.5,sigma)-(Y-y0+0.5)*self.gauss(Y,y0-0.5,sigma))/(sigma**3)).flatten()
        nk = H.flatten()
        muk = mu.flatten()
        Dl = np.divide(nk,muk,out=np.zeros_like(nk,dtype=np.float32),where=muk!=0)
        Dl2 = np.divide(nk,muk**2,out=np.zeros_like(nk,dtype=np.float32),where=muk!=0)
        Dl2[Dl2 == np.inf] = 0
        Hl = np.zeros((2,2))
        Hl[0,0] = np.sum(Hm[:,0,0]*Dl) + np.sum(Dl2*Jm[:,0]*Jm[:,0])
        Hl[0,1] = np.sum(Hm[:,0,1]*Dl) + np.sum(Dl2*Jm[:,0]*Jm[:,1])
        Hl[1,0] = np.sum(Hm[:,1,0]*Dl) + np.sum(Dl2*Jm[:,1]*Jm[:,0])
        Hl[1,1] = np.sum(Hm[:,1,1]*Dl) + np.sum(Dl2*Jm[:,1]*Jm[:,1])
        return Hl
        
    def gauss(self,x,x0,sigma):
        return np.exp(-((x-x0)**2)/(2*sigma**2))


class NoisyGaussianPSFMLE:
    def __init__(self,X):
        self.X = X
    def estimate(self,theta0,alpha,gain,I0,eta,rmu,rvar):
        L = self.X.shape[0]
        def neg_loglike(theta):
            x0,y0 = theta
            x = np.arange(0,L)
            y = np.arange(0,L)
            X,Y = np.meshgrid(x,y)
            intx = 0.5*(erf((X+0.5-x0)*alpha)-erf((X-0.5-x0)*alpha))
            inty = 0.5*(erf((Y+0.5-y0)*alpha)-erf((Y-0.5-y0)*alpha))
            pmu = gain*I0*eta*intx*inty
            mu = rmu + pmu
            var = rvar + pmu
            nll = 0.5*np.log(2*np.pi*var) + ((self.X-mu)**2)/(2*var)
            ll = np.sum(nll)
            return ll
        jacobian_ = jacobian(neg_loglike)
        bounds = [(4,6),(4,6)]
        res1 = minimize(neg_loglike, theta0, method = 'BFGS', \
        	       options={'disp': True}, bounds=bounds, jac = jacobian_)
        return res1

    def plot_mle_surf(self,alpha,gain,I0,eta,rmu,rvar,nx=100,ny=100):
        L = self.X.shape[0]
        def neg_loglike(theta):
            x0,y0 = theta
            x = np.arange(0,L)
            y = np.arange(0,L)
            X,Y = np.meshgrid(x,y)
            intx = 0.5*(erf((X+0.5-x0)*alpha)-erf((X-0.5-x0)*alpha))
            inty = 0.5*(erf((Y+0.5-y0)*alpha)-erf((Y-0.5-y0)*alpha))
            pmu = gain*I0*eta*intx*inty
            mu = rmu + pmu
            var = rvar + pmu
            nll = 0.5*np.log(2*np.pi*var) + ((self.X-mu)**2)/(2*var)
            ll = np.sum(nll)
            return ll
        xs = np.linspace(0,L,nx); ys = np.linspace(0,L,ny)
        xv, yv = np.meshgrid(xs,ys)
        lsurf = np.zeros((nx,ny))
        jtensor = np.zeros((nx,ny,2))
        htensor = np.zeros((nx,ny,2,2))
        jacobian_ = jacobian(neg_loglike)
        hessian_ = hessian(neg_loglike)
        for i in range(nx):
            for j in range(ny):
                theta = np.append(xv[i,j],yv[i,j])
                lsurf[i,j] = neg_loglike(theta)
                jtensor[i,j] = jacobian_(theta)
                htensor[i,j] = hessian_(theta)
        jsurf = np.sqrt(np.sum(jtensor**2,axis=-1))
        hsurf = np.trace(htensor,axis1=2,axis2=3)
        fig, ax = plt.subplots(2,2,figsize=(4,3))
        im1 = ax[0,0].imshow(self.X,cmap='plasma')
        ax[0,0].set_xticks([]);ax[0,0].set_yticks([])
        plt.colorbar(im1, ax=ax[0,0], label='ADU')
        im2 = ax[0,1].imshow(lsurf,cmap='plasma')
        ax[0,1].set_xticks([]);ax[0,1].set_yticks([])
        plt.colorbar(im2, ax=ax[0,1], label=r'$\ell(x,y)$')
        im3 = ax[1,0].imshow(jsurf,cmap='plasma')
        ax[1,0].set_xticks([]);ax[1,0].set_yticks([])
        plt.colorbar(im3, ax=ax[1,0], label=r'$|\nabla \ell|$')
        im4 = ax[1,1].imshow(hsurf,cmap='plasma')
        ax[1,1].set_xticks([]);ax[1,1].set_yticks([])
        plt.colorbar(im4, ax=ax[1,1], label=r'$\mathrm{Tr}(H)$')
        plt.tight_layout()
        plt.show()

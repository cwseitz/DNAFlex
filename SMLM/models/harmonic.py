from __future__ import division
import numpy as np
from scipy.linalg import logm

class HarmonicOscillator:
    def __init__(self, N, m, T, D, k, kB, fs, dt, g, w0, w, tau ):
        self.N   =  N           # # of data points
        self.m   =  m           # mass
        self.T   =  T           # temperature
        self.D   =  D           # diffusion constant
        self.k   =  k           # stiffness
        self.kB  =  kB          # Boltzmann constant
        self.fs  =  fs          # sampling frequency
        self.dt  =  dt          # step size
        self.g   =  g           # gamma
        self.w0  =  w0          # natural frequency
        self.w   =  w           # frequency
        self.tau =  tau         # relaxation time

        
    def calcSigma(self):
        kmw  = self.k/(self.m*self.w*self.w);   
        dtbt = -self.dt/self.tau; ee=np.exp(dtbt);  
        dd   = self.D/(self.w*self.w*self.tau*self.tau);  
        tt   = self.w*self.dt;   cc  = np.cos(tt);  ss=np.sin(tt)
                
        s1 = (self.kB*self.T/self.k)*(1-ee*(kmw*ss*ss+(cc+ss/(2*self.w*self.tau))**2))
        s2 = dd*ee*ss*ss
        s3 = (self.kB*self.T/self.m)*(1-ee*(kmw*ss*ss+(cc-ss/(2*self.w*self.tau))**2)) 
        return s1, s3, s2


    def calcLambda(self):
        ii  = np.eye(2)
        ll = np.asanyarray([[0, -1], [self.k/self.m, self.g/self.m]])
        ee = np.exp(-self.dt/(2*self.tau)); wt2=2*self.w*self.tau
        cc=np.cos(self.w*self.dt); ss=np.sin(self.w*self.dt) 

        Lambda = ee*((cc+ss/wt2)*ii - ll*ss/self.w ) 
        return np.real(Lambda)


    def calcXV(self, Lambda, Sigma):
        x = np.zeros([self.N,1])
        v = np.zeros([self.N,1])

        s1, s3, s2 = Sigma

        for j in np.arange(0,self.N-1):
            oldvec = np.array([x[j],v[j]])
            randgauss = np.random.randn(2,1)
            delx = np.sqrt(s1)*randgauss[0]
            delv = (s2/(np.sqrt(s1)))*randgauss[0]+(np.sqrt(s3 - ((s2**2)/(s1))))*randgauss[1]
            delvec = np.array([delx,delv])
            updatevec = np.dot(Lambda,oldvec)+delvec
            x[j+1] = updatevec[0]
            v[j+1] = updatevec[1]
        return x,v


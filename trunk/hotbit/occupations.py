import numpy as nu
from scipy.optimize import brentq

class Occupations:
    def __init__(self,nel,width):
        self.nel=nel
        self.width=width
        
    def fermi(self,e,mu,width):
        y=(e-mu)/width
        if y>100:
            return 0.0
        elif y<-100:
            return 2.0
        else:
            return 2.0/(nu.exp(y)+1)        
                
    def root_function(self,mu):
        args=(self.e-mu)/self.width
        exps=nu.exp(args)
        #exps=where(args>100,exps,0.0)
        #exps=where(args<-100,2.0,exps)
        occus=2/(exps+1)                
        return occus.sum()-self.nel
        
    def occupy(self,e,width=None):
        if width!=None:
            self.width=width
        self.e=e
        if nu.mod(self.nel,2)==0:
            mu_guess=self.e[int(round(self.nel/2.0))]
        else:
            mu_guess=0.5*( self.e[int(nu.floor(self.nel/2.0))]+self.e[int(nu.ceil(self.nel/2.0))] )
        dmu=self.width
        try:
            self.mu=brentq(self.root_function,mu_guess-dmu,mu_guess+dmu)
        except:
            self.mu=brentq(self.root_function,self.e[0]-dmu,self.e[-1]+dmu)
        self.f=nu.array([self.fermi(ek,self.mu,self.width) for ek in self.e])
        return self.f
        
    def plot(self):
        import pylab as pl
        pl.plot(self.e,self.f)
        pl.scatter(self.e,self.f)
        pl.title('occupations')
        pl.ylim(-0.2,2.2)
        pl.xlabel('energy (Ha)')
        pl.ylabel('occupation')
        pl.show()
        
            
        
        
    
        
       
from box.interpolation import SplineFunction
import numpy as np

class ArbitraryDistribution:
    """ Class for generating random numbers form given distribution. """
    
    def __init__(self,distr=None,xgrid=None):
        """ 
        Initialize either with given probability distribution 'distr(x)',
        or give the distribution on grid (given grid or 0,1,2,... by default)
        """
        try:
            distr(0.0)
            self.f=distr
        except:
            self.x=xgrid
            if self.x is None:
                self.x=range(len(distr))
            self.distr=distr
            self.f=SplineFunction(self.x,self.distr,k=1)
                
        # put the probability distribution f(x) on grid
        self.N=1000
        self.lim=self.f.limits()
        self.xgrid=np.linspace(self.lim[0],self.lim[1],self.N)
        self.fvals=[self.f(z) for z in self.xgrid]
        
        # construct the distribution function F(x)
        cum=[0.0]
        for i in range(self.N-1):
            cum.append(cum[-1]+self.fvals[i]+self.fvals[i+1])
        cum=np.array(cum)/2*(self.lim[1]-self.lim[0])/self.N
        self.norm=cum[-1]
        cum=cum/self.norm
        self.F=SplineFunction(self.xgrid,cum,k=1)
        
    def __call__(self):
        """ Generate random number from given distribution. """
        rnd=np.random.uniform(0.0,1.0)
        return self.F.solve(rnd)
    
    def verify(self,nr=10000):
        """ Plot and compare given distribution and simulated one. """
        c=np.array([self() for rnd in range(nr)])
        
        import pylab as pl
        pl.hist(c,len(self.xgrid),normed=True)
        distr=[self.f(x)/self.norm for x in self.xgrid]
        pl.plot(self.xgrid,distr,c='r')
        pl.show()
            
if __name__=='__main__':
    y=[0,1,2,3,0.1,0.1,0.1,5]
    f=ArbitraryDistribution(distr=y)
    f.verify()        
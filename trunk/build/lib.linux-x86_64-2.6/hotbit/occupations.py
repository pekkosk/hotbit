import numpy as np
from scipy.optimize import brentq

class Occupations:
    def __init__(self,nel,width,wk):
        '''
        Initialize parameters for occupations.
        
        @param nel: Number of electrons 
        @param width: Fermi-broadening
        @param wk: k-point weights
        '''
        self.nel=nel
        self.width=width
        self.wk=wk
        self.nk=len(wk)
        
        
    def get_mu(self):
        """ Return the Fermi-level (or chemical potential). """
        return self.mu


    def fermi(self,mu):
        """ Occupy states with given chemical potential. 
        
        Occupations are 0...2; without k-point weights
        """
        args = (self.e-mu)/self.width
        exps = np.exp(args)
        return 2/(exps+1)        
        

    def root_function(self,mu):
        """ This function is exactly zero when mu is right. """
        # w_k transpose(f_ka) = w_k f_ak -> transpose -> f_ka*w_k
        f = self.fermi(mu)
        kf = (self.wk * f.transpose()).transpose()
        return kf.sum()-self.nel


    def occupy(self,e):
        '''
        Calculate occupation numbers with given Fermi-broadening.
        
        @param e: e[k,a] energy of k-point, state a
        @param wk: wk[:] weights for k-points
        @param width: The Fermi-broadening
        '''
        self.e=e
        ef = np.sort( e.flatten() )
        try:
            # make the first guess (simple guess, assuming equal k-point weights)
            
            guess = self.ef[int(round(self.nk*self.nel/2.0))]
            dmu = self.width*5
            mu  = brentq(self.root_function,guess-dmu,guess+dmu)
        except:
            # probably a bad guess
            dmu = self.width 
            mu = brentq(self.root_function,ef[0]-dmu,ef[-1]+dmu)    
        f = self.fermi(mu)
        self.mu, self.f = mu, f
        return f


    def plot(self):
        import pylab as pl
        for ik in range(self.nk):
            pl.plot(self.e[ik,:],self.f[ik,:])
            pl.scatter(self.e[ik,:],self.f[ik,:])
        pl.title('occupations')
        pl.xlabel('energy (Ha)')
        pl.ylabel('occupation')
        pl.show()








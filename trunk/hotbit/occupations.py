import numpy as np
from scipy.optimize import brentq
import sys
if sys.version_info < (2,6):
    MAX_EXP_ARGUMENT = np.log(1E90)
else:
    MAX_EXP_ARGUMENT = np.log(sys.float_info.max)


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
        # np.exp will return np.inf with RuntimeWarning if the input value
        # is too large, better to feed it nu.inf in the beginning
        args = np.where(args < MAX_EXP_ARGUMENT, args, np.inf)
        exps = np.exp(args)
        return 2/(exps+1)        
        
        
    def fermi_function(self,e):
        """
        Return Fermi function for given energy [0,1]
        """
        arg = (e-self.mu)/self.width
        return 1/(np.exp(arg)+1)  
        

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
        eflat = e.flatten()
        ind = np.argsort( eflat )
        e_sorted = eflat[ind]
        n_sorted = (self.wk*np.ones_like(e).transpose()*2).transpose().flatten()[ind]
        sum = n_sorted.cumsum() 
        ifermi = np.searchsorted(sum,self.nel)
        
        xtol = 1E-13
        
        try:
            # make the first guess (simple guess, assuming equal k-point weights)
            elo = e_sorted[ifermi-1]
            ehi = e_sorted[ifermi-1]
            guess = e_sorted(ifermi)
            dmu = np.max((self.width,guess-elo,ehi-guess))
            mu  = brentq(self.root_function,guess-dmu,guess+dmu,xtol=xtol)
        except:
            # probably a bad guess        
            dmu = self.width 
            mu = brentq(self.root_function,e_sorted[0]-dmu,e_sorted[-1]+dmu,xtol=xtol)
        
        if np.abs( self.root_function(mu) )>1E-10:
            raise RuntimeError('Fermi level could not be assigned reliably. Has the system fragmented?')
            
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








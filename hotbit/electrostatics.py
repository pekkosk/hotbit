# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

import numpy as nu
from scipy.special import erf
from hotbit import auxil
from math import sqrt,pi,exp
from weakref import proxy
dot=nu.dot


class Electrostatics:
    def __init__(self,calc):
        self.calc=proxy(calc)
        self.dq=nu.zeros(len(calc.el))
        self.norb=calc.el.get_nr_orbitals()
        self.SCC=calc.get('SCC')
        self.N=len(calc.el)
        self.timer=calc.timer
        self.cut=calc.ia.get_cut()
        self.gamma_potential=None
        
    def __del__(self):
        print "Electrostatics deleted"

    def __call__(self):
        """ Electrostatic (TB) potential at r. """
        return None
    
    def set_dq(self,dq):
        """ (Re)setting dq gives new gamma-potential. """
        if nu.all(abs(dq-self.dq)<1E-13) and self.gamma_potential!=None:
            return
        self.dq=dq
        self.gamma_potential=nu.zeros((self.N,))
        self.gamma_potential_der=nu.zeros((self.N,3))
        for i in range(self.N):
            self.gamma_potential[i]=dot(self.gamma_table[i,:],-self.dq[:])
            for a in range(3):
                self.gamma_potential_der[i,a]=dot(self.gamma_der[i,:,a],-self.dq[:])      
        
    def coulomb_energy(self):
        """ Return Coulomb energy with given Mulliken charges. """
        if not self.SCC:
            return 0.0
        else:
            return dot(-self.dq[:],self.gamma_potential[:])/2
        
    def gamma_forces(self):
        """ Return forces due to electrostatic interactions. """
        
        if not self.SCC:
            return nu.zeros((self.N,3))
        else:
            self.timer.start('f_es')
            f=[]
            for i in range(self.N):
                f.append( (-self.dq[i])*self.gamma_potential_der[i] )                                                
            self.timer.stop('f_es')                
            return nu.array(f)
            
        
    def construct_H1(self,dq=None):
        """ Make the electrostatic part of the Hamiltonian. """
        self.timer.start('electrostatic H')
        if dq!=None:
            self.set_dq(dq)        
        H1=nu.zeros((self.norb,self.norb))        
        
        gp=self.gamma_potential
        ext=self.ext
        for i,j,o1i,o1j,noi,noj in self.pairs:
            # internal electrostatics
            phi_ij=(-1)*(gp[i]+gp[j])/2
            H1[o1i:o1i+noi,o1j:o1j+noj]=phi_ij
            H1[o1j:o1j+noj,o1i:o1i+noi]=phi_ij
            # external electrostatics
            ext_ij=(-1)*(ext[i]+ext[j])/2
            H1[o1i:o1i+noi,o1j:o1j+noj]+=ext_ij
            H1[o1j:o1j+noj,o1i:o1i+noi]+=ext_ij    
    
        self.timer.stop('electrostatic H')                        
        self.H1=H1                               
        return self.H1
                    
    def get_H1(self):
        """ Get the current electrostatic Hamiltonian. """
        return self.H1    
                    
    def construct_tables(self):
        """ Stuff calculated _once_ for given coordinates. """
        self.timer.start('es tables')
        el = self.calc.el
        #7: return (i,o1i,noi,j,o1j,noj) 
        self.pairs=el.get_ia_atom_pairs(['i','j','o1i','o1j','noi','noj']) # this also only once
        #print [self.el.orbitals(i,indices=True)[0] for i in range(self.N)]
        #print [self.el.orbitals(i,number=True) for i in range(self.N)]
        #self.orb1=nu.array([self.el.orbitals(i,indices=True)[0] for i in range(self.N)])
        #self.norbs=nu.array([self.el.orbitals(i,number=True) for i in range(self.N)] )
        
        g=nu.zeros((self.N,self.N))
        dg=nu.zeros((self.N,self.N,3))
        for i,j in el.get_ia_atom_pairs(['i','j']):
            g[i,j]=self.gamma(i,j)
            g[j,i]=g[i,j]
            if i!=j:
                dg[i,j,:]=self.gamma(i,j,der=1)
                dg[j,i,:]=-dg[i,j,:]
        self.gamma_table, self.gamma_der=g, dg            
        self.ext = [self.calc.env.phi(i) for i in range(self.N)]
        self.timer.stop('es tables')
            
    def gamma(self,i,j,der=0):
        """ Return the gamma function for atoms i and j. der=1 for gradient. """
        el = self.calc.el
        r=el.distance(i,j)
        cut=self.calc.get('gamma_cut')
        ei, ej=[el.get_element(k) for k in [i,j]]
        wi, wj=ei.get_FWHM(), ej.get_FWHM()
        const=2*nu.sqrt( nu.log(2.0)/(wi**2+wj**2) )
        if i==j:
            assert der==0
            return ei.get_U()
        else:
            ecr=erf(const*r) 
            if der==0:
                if cut==None:
                    return ecr/r
                else:
                    return ecr/r*(1-erf(r/cut))
            elif der==1:
                vec=el.vector(i,j)
                decr=2/sqrt(pi)*exp(-(const*r)**2)*const
                if cut==None:
                    return (decr - ecr/r)*vec/r**2         
                else:          
                    f=(1-erf(r/cut))    
                    df=-2/sqrt(pi)*exp(-(r/cut)**2)/cut
                    return (df*ecr + f*decr - f*ecr/r)*vec/r**2  

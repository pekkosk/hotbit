# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

import numpy as nu
from scipy.special import erf
from hotbit import auxil
from math import sqrt,pi,exp
from weakref import proxy
from math import sqrt
norm=nu.linalg.norm
dot=nu.dot


class Electrostatics:
    def __init__(self,calc):
        self.calc=proxy(calc)
        self.norb=calc.el.get_nr_orbitals()
        self.SCC=calc.get('SCC')
        self.N=len(calc.el)
        self.dq=nu.zeros((self.N))
        self.G=nu.zeros((self.N,self.N))
        self.dG=nu.zeros((self.N,self.N,3))


    def set_dq(self,dq):
        """ (Re)setting dq gives new gamma-potential. """
        self.dq=dq
        self.epsilon = nu.array( [dot(self.G[i,:],self.dq) for i in range(self.N)] )
            

    def coulomb_energy(self):
        """ Return Coulomb energy. """
        self.calc.start_timing('ecoul')
        ecoul = 0.0
        if not self.SCC:
            ecoul = 0.0
        else:
            ecoul = 0.5 * dot(self.dq,self.epsilon)
        self.calc.stop_timing('ecoul')
        return ecoul


    def gamma_forces(self):
        """ Return forces due to electrostatic interactions. """
        if not self.SCC:
            return nu.zeros((self.N,3))
        else:
            self.calc.start_timing('f_es')
            f = nu.zeros((self.N,3))
            for a in range(3):
                depsilon = [dot(self.dq,self.dG[k,:,a]) for k in range(self.N)]
                f[:,a] = self.dq * depsilon 
            self.calc.stop_timing('f_es')
            return f


    def construct_h1(self,dq=None):
        """ Make the electrostatic part of the Hamiltonian. """
        self.calc.start_timing('h1')
        if dq!=None:
            self.set_dq(dq)
        h1=nu.zeros((self.norb,self.norb))

        ext=self.ext            
        lst = self.calc.el.get_property_lists(['i','o1','no'])
        for i,o1i,noi in lst:
            for j,o1j,noj in lst:
                # internal electrostatics from Mulliken charges
                eij = 0.5 * (self.epsilon[i] + self.epsilon[j])
                h1[o1i:o1i+noi,o1j:o1j+noj] = eij
                h1[o1j:o1j+noj,o1i:o1i+noi] = eij
                # external electrostatics 
                # NOTE: ext is the electrostatic potential, so we have minus-sign!
                ext_ij = 0.5 * (-1)*(ext[i]+ext[j])
                h1[o1i:o1i+noi,o1j:o1j+noj] += ext_ij
                h1[o1j:o1j+noj,o1i:o1i+noi] += ext_ij
                
        self.calc.stop_timing('h1')
        self.h1=h1
        return self.h1


    def get_h1(self):
        """ Get the current electrostatic Hamiltonian. """
        return self.h1


    def construct_Gamma_matrix(self):
        '''
        Construct the G-matrix and its derivative.
        
        Done once for each geometry; 
        G_ij = sum_n gamma_ij(Rijn)
        dG_ij = sum_n gamma'_ij(Rijn) hat(Rijn) 
        '''
        self.calc.start_timing('gamma matrix')
        g=self.gamma
        G=nu.zeros((self.N,self.N))
        dG=nu.zeros((self.N,self.N,3))
        
        lst=self.calc.el.get_property_lists(['i','s'])
        for i,si in lst:
            for j,sj in lst:
                Rijn = self.calc.el.Rn[:,j,:] - self.calc.el.Rn[0,i,:]
                for n,rij in enumerate(Rijn):
                    dij = sqrt( rij[0]**2+rij[1]**2+rij[2]**2 )
                    G[i,j] = G[i,j] + g(si,sj,dij)
                    if i==j and n==0: continue
                    dG[i,j,:] = dG[i,j,:] + g(si,sj,dij,der=1)*rij/dij 
        self.G, self.dG = G, dG
                
        self.ext = [self.calc.env.phi(i) for i in range(self.N)]
        self.calc.stop_timing('gamma matrix')

                            


    def gamma(self,si,sj,r,der=0):
        '''
        Return the gamma for atoms i and j.
        
        gamma_ij(r) = erf(c*r)/r * erfc(r/cut)
        
        @param i: atom index (only element information)
        @param j: atom index2 (only element information)
        @param r: distance between atoms.
        @param der: 0 = gamma(r); 1 = dgamma(r)/dr 
        '''
        cut = self.calc.get('gamma_cut') # for more rapid decay of Coulomb potential
        ei, ej = [self.calc.el.get_element(s) for s in (si,sj)]
        wi, wj = ei.get_FWHM(), ej.get_FWHM()
        const=2*nu.sqrt( nu.log(2.0)/(wi**2+wj**2) )
        if r<1E-10:
            assert (der==0 and si==sj)
            return ei.get_U()
        else:
            ecr=erf(const*r)
            if der==0:
                if cut==None:
                    return ecr/r
                else:
                    return ecr/r*(1-erf(r/cut))
            elif der==1:
                decr=2/sqrt(pi)*exp(-(const*r)**2)*const
                if cut==None:
                    return (decr - ecr/r)/r 
                else:
                    f=(1-erf(r/cut))
                    df=-2/sqrt(pi)*exp(-(r/cut)**2)/cut
                    return (df*ecr + f*decr - f*ecr/r)/r 
            
                

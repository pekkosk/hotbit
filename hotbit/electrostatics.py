"""
Compute electrostatic interaction and construct the 'electrostatic'
matrix h1 containing the shift in on-site energies.

This module only adds the contribution of the shape of the charge (currently
only Gaussian). The interaction of point charges is compute by either the Ewald
or MultipoleExpansion class.
"""

# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

from math import exp, log, pi, sqrt
from weakref import proxy
from math import sqrt

import numpy as nu
from scipy.special import erf, erfc
norm=nu.linalg.norm
dot=nu.dot

from ase.units import Bohr

from hotbit.ewald_sum import EwaldSum
from hotbit.multipole_expansion import MultipoleExpansion


# Some constants
log2 = log(2.0)


# List of available point-charge Coulomb solvers
str_to_solver = {
    'ewald': EwaldSum,
    'me':    MultipoleExpansion
    }


class Electrostatics:
    def __init__(self, calc, solver=None, accuracy_goal=12, solver_args={}):
        self.calc=proxy(calc)
        self.norb=calc.el.get_nr_orbitals()
        self.SCC=calc.get('SCC')
        self.N=len(calc.el)
        self.dq=nu.zeros((self.N))
        self.G=nu.zeros((self.N,self.N))
        self.dG=nu.zeros((self.N,self.N,3))

        self.accuracy_goal = accuracy_goal

        if solver is None:
            self.solver = None
            self.gamma  = self.gamma_direct
        else:
            self.solver = str_to_solver[solver](timer = self.calc.timer,
                                                **solver_args)
            self.gamma  = self.gamma_correction


    def set_dq(self,dq):
        """ (Re)setting dq gives new gamma-potential. """
        self.dq=dq
        self.epsilon = nu.sum(self.G*self.dq, axis=1)
        if self.solver is not None:
            self.solver.update(self.calc.el.atoms, dq)
            self.epsilon += self.solver.get_potential()*Bohr


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
            if self.solver is None:
                E = nu.zeros((self.N,3))
            else:
                E = self.solver.get_field()
            depsilon  = nu.sum(self.dq.reshape(1, -1, 1)*self.dG, axis=1)
            f         = self.dq.reshape(-1, 1) * ( depsilon + E )
            self.calc.stop_timing('f_es')
            return f


    def construct_h1(self,dq=None):
        """ Make the electrostatic part of the Hamiltonian. """
        self.calc.start_timing('h1')
        if dq!=None:
            self.set_dq(dq)
                          
        lst = self.calc.el.get_property_lists(['i','o1','no'])
        
        aux = nu.zeros((self.norb,self.norb))
        for i,o1i,noi in lst:
            aux[o1i:o1i+noi,:] = self.epsilon[i]        
        h1 = 0.5 * ( aux+aux.transpose() )
                
        # external electrostatics
        if nu.any(abs(self.ext)>1E-12):
            aux = nu.zeros((self.norb,self.norb))
            for i,o1i,noi in lst:
                aux[o1i:o1i+noi,:] = self.ext[i]
            h1 = h1 + 0.5 * (-1) * (aux+aux.transpose())
                
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

        # Heuristics to estimate the cut-off for the Coulomb correction.
        # Not yet used.
        if self.gamma == self.gamma_direct:
            cutoff  = self.calc.get('gamma_cut')
        else:
            min_U = min([ self.calc.el.get_element(s).get_U()
                          for s in self.calc.el.get_present() ])
            cutoff = sqrt(log(10.0)*self.accuracy_goal/(sqrt(pi/2)*min_U))

        g=self.gamma
        G=nu.zeros((self.N,self.N))
        dG=nu.zeros((self.N,self.N,3))
        rijn, dijn = self.calc.el.get_distances()

        lst=self.calc.el.get_property_lists(['i','s'])
        
        for i,si in lst:
            for j,sj in lst[i:]:
                G[i,j] = sum( [g(si,sj,d) for d in dijn[:,i,j]] )

                aux = nu.array([ g(si,sj,d,der=1)*r/d
                                 for (d,r) in zip(dijn[:,i,j],rijn[:,i,j]) ])
                if i==j:
                    # exclude n=(0,0,0) from derivatives
                    dG[i,j] = aux[1:,:].sum(axis=0)
                elif i!=j:
                    G[j,i] = G[i,j] 
                    dG[i,j] = aux.sum(axis=0)
                    dG[j,i] = -dG[i,j]

        self.G, self.dG = G, dG
                
        self.ext = nu.array( [self.calc.env.phi(i) for i in range(self.N)] )
        self.calc.stop_timing('gamma matrix')


    def gamma_correction(self,si,sj,r,der=0):
        '''
        Return the gamma for atoms i and j.
        Note that this does not contain the contribution of
        the point charges, and is hence short ranged.

        gamma_ij(r) = -erfc(c*r)/r * erfc(r/cut)
        
        @param i: atom index (only element information)
        @param j: atom index2 (only element information)
        @param r: distance between atoms.
        @param der: 0 = gamma(r); 1 = dgamma(r)/dr 
        '''
        # If the correction is used, the system is periodic anyway and hence
        # 'gamma_cut' is not supported
        ei, ej = [self.calc.el.get_element(s) for s in (si,sj)]
        wi, wj = ei.get_FWHM(), ej.get_FWHM()
        const=2*sqrt( log2/(wi**2+wj**2) )
        if r<1E-10:
            assert si==sj
            if der==1:
                return nu.zeros((3))
            else:
                return ei.get_U()
        else:
            ecr = -erfc(const*r)
            if der==0:
                return ecr/r
            elif der==1:
                decr=2/sqrt(pi)*exp(-(const*r)**2)*const - 1/r
                return (decr - ecr/r)/r 
            
                
    def gamma_direct(self,si,sj,r,der=0):
        '''
        Return the gamma for atoms i and j.

        gamma_ij(r) = -erfc(c*r)/r * erfc(r/cut)
        
        @param i: atom index (only element information)
        @param j: atom index2 (only element information)
        @param r: distance between atoms.
        @param der: 0 = gamma(r); 1 = dgamma(r)/dr 
        '''
        # for more rapid decay of Coulomb potential
        cut = self.calc.get('gamma_cut')
        ei, ej = [self.calc.el.get_element(s) for s in (si,sj)]
        wi, wj = ei.get_FWHM(), ej.get_FWHM()
        const=2*sqrt( log2/(wi**2+wj**2) )
        if r<1E-10:
            assert si==sj
            if der==1:
                return nu.zeros((3))
            else:
                return ei.get_U()
        else:
            ecr = erf(const*r)
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
                    f=erfc(r/cut)
                    df=-2/sqrt(pi)*exp(-(r/cut)**2)/cut
                    return (df*ecr + f*decr - f*ecr/r)/r 
            
                

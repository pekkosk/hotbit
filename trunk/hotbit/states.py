# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

from solver import Solver
from hotbit.fortran.misc import fortran_rho
from hotbit.fortran.misc import fortran_rho0
from hotbit.fortran.misc import fortran_rhoe0
from hotbit.fortran.misc import symmetric_matmul
from hotbit.fortran.misc import matmul_diagonal
from hotbit.fortran.misc import fortran_fbs
from occupations import Occupations
import numpy as nu

class States:
    
    def __init__(self,calc,timer,elements,interactions):
        self.calc=calc
        self.timer=timer
        self.el=elements
        self.ia=interactions
        self.solver=Solver(calc,timer)
        width=calc.get('width')
        self.occu=Occupations(self.el.get_number_of_electrons(),width=width)
        self.nat=len(elements)
        self.es=calc.es
        self.norb=self.el.get_nr_orbitals()
        self.prev_dq=[None,None]
        self.count=0
        self.SCC=calc.get('SCC')
        self.rho=None
        self.rhoe=None
        
    def guess_dq(self):
        n=len(self.el)
        if not self.SCC:
            return nu.zeros((n,))
        if self.count==0:
            return nu.zeros((n,)) - self.calc.get('charge')/n
        elif self.count==1:    # use previous charges 
            return self.prev_dq[0]
        else:                  # linear extrapolation
            return self.prev_dq[0] + (self.prev_dq[0]-self.prev_dq[1])
        
    def solve(self):
        self.timer.start('solve')
        dq=self.guess_dq()
        self.norb=self.el.get_nr_orbitals()
        self.H0, self.S, self.dH0, self.dS=self.ia.get_matrices()
        self.e, self.wf=self.solver.get_states(self,dq,self.H0,self.S)
        self.el.set_solved()
        self.large_update()
        self.count+=1
        self.timer.stop('solve')
    
    def update(self,e,wf):
        """ Update all essential stuff from given wave functions. """
        self.timer.start('update')
        self.e=e
        self.wf=wf
        self.f=self.occu.occupy(e)
        self.timer.start('rho')
        self.rho0=fortran_rho0(self.wf,self.f,self.el.nr_ia_orbitals,self.el.ia_orbitals,self.norb)
        self.timer.stop('rho')
        self.rho0S_diagonal=matmul_diagonal(self.rho0,self.S,self.norb)
        if self.SCC:
            self.dq=self.mulliken()
            self.es.set_dq(self.dq)
        self.timer.stop('update')       
        
    def large_update(self):
        """ Update stuff from eigenstates needed later for forces etc. """
        self.timer.start('final update')
        self.Swf=None
        #self.rho=fortran_rho(self.wf,self.f,self.norb) # the complete density matrix (even needed?)
        self.dH=nu.zeros_like(self.dH0)
        if self.SCC:
            self.prev_dq=[self.dq, self.prev_dq[0]]
            for a in range(3):
                self.dH[:,:,a]=self.dH0[:,:,a] + self.es.get_H1()*self.dS[:,:,a]
        else:
            self.dH=self.dH0                 
            
        # density matrix weighted by eigenenergies 
        self.rhoe0=fortran_rhoe0(self.wf,self.f,self.e,self.el.nr_ia_orbitals,self.el.ia_orbitals,self.norb)                     
        self.timer.stop('final update')
                                   
    def get_dq(self):
        return self.dq
            
    def get_eigenvalues(self):
        return self.e
    
    def get_occupations(self):
        return self.f
    
    def get_homo(self):
        """ Return highest orbital with occu>0.99. """
        for i in range(self.norb)[::-1]:
            if self.f[i]>0.99: return i
        
    def get_lumo(self):
        """ Return lowest orbital with occu<1.01. """
        for i in range(self.norb):
            if self.f[i]<1.01: return i
        
    def get_hoc(self):
        """ Return highest partially occupied orbital. """
        for i in range(self.norb)[::-1]:
            if self.f[i]>1E-9: return i
                        
    def mulliken(self):
        """ Return excess Mulliken populations. """
        q=[]
        for i in range(self.nat):
            orbitals=self.el.orbitals(i,indices=True)
            q.append( sum(self.rho0S_diagonal[orbitals]) )           
        return nu.array(q)-self.el.get_valences()
            
    def mulliken_transfer(self,k,l):
        """ Return Mulliken transfer charges between states k and l. """
        if self.Swf==None:
            self.Swf=symmetric_matmul(self.S,self.wf)
        q=[]
        for i in range(self.nat):            
            iorb=self.el.orbitals(i,indices=True)
            qi=sum( [self.wf[oi,k]*self.Swf[oi,l]+self.wf[oi,l]*self.Swf[oi,k] for oi in iorb] )
            q.append(qi/2)
        return nu.array(q)
                
    def band_structure_energy(self):
        """ Return band structure energy. """
        return nu.trace(nu.dot(self.rho0,self.H0)) 
                           
    def band_structure_forces(self):
        """ Return forces arising from band structure. """
        self.timer.start('f_bs')
        
        norbs=self.el.nr_orbitals
        inds=self.el.atom_orb_indices2
        
        f=fortran_fbs(self.rho0,self.rhoe0,self.dH,self.dS,norbs,inds,self.norb,self.nat)
        self.timer.stop('f_bs')                
        return -2*nu.real(f)
                                        
      
 

# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

from solver import Solver
from electrostatics import Electrostatics
from occupations import Occupations
from hotbit.fortran.misc import fortran_rho
from hotbit.fortran.misc import fortran_rho0
from hotbit.fortran.misc import fortran_rhoe0
from hotbit.fortran.misc import fortran_rhoc
from hotbit.fortran.misc import fortran_rhoec
from hotbit.fortran.misc import symmetric_matmul
from hotbit.fortran.misc import matmul_diagonal
from hotbit.fortran.misc import fortran_fbs
from hotbit.fortran.misc import fortran_fbsc
import numpy as nu
pi=nu.pi
from weakref import proxy

class States:

    def __init__(self,calc):
        self.es=Electrostatics(calc)
        self.solver=Solver(calc)
        self.calc=proxy(calc)
        width=calc.get('width')
        
        self.nat=len(calc.el)
        self.norb=calc.el.get_nr_orbitals()
        self.prev_dq=[None,None]
        self.count=0
        self.SCC=calc.get('SCC')
        self.rho=None
        self.rhoe0=None
        
        kpts=calc.get('kpts')
        if isinstance(kpts,tuple):
            # set up equal-weighted and spaced k-point mesh
            if 0 in kpts:
                raise AssertionError('Each direction must have at least one k-point! (Gamma-point)')
            
            self.kl=[]
            for i in range(3):
                spacing = 2*pi/kpts[i]
                self.kl.append( nu.linspace(-pi+spacing/2,pi-spacing/2,kpts[i]) )
            
            self.k=[]    
            for a in range(kpts[0]):
                for b in range(kpts[1]):
                    for c in range(kpts[2]):
                        self.k.append( nu.array([self.kl[0][a],self.kl[1][b],self.kl[2][c]]) )
            self.nk=nu.prod(kpts)
            self.wk=nu.ones(self.nk)/self.nk
            
        else:
            # work with a given set of k-points
            self.nk=len(kpts)
            self.k=kpts
            self.wk=nu.ones(self.nk)/self.nk
            self.kl=None
            
        assert sum(self.wk)-1.0<1E-13
        pbc=self.calc.el.pbc
        for i in range(3):
            for k in self.k:
                if k[i]>1E-10 and not pbc[i]:
                    raise AssertionError('Do not set (non-zero) k-points in non-periodic direction!')
            
        self.occu=Occupations(calc.el.get_number_of_electrons(),width,self.wk)
        # TODO: check what k-weights add up to one

    def __del__(self):
        pass


    def guess_dq(self):
        n=len(self.calc.el)
        if not self.SCC:
            return nu.zeros((n,))
        if self.count==0:
            return nu.zeros((n,)) - float(self.calc.get('charge'))/n
        elif self.count==1:    # use previous charges
            return self.prev_dq[0]
        else:                  # linear extrapolation
            return self.prev_dq[0] + (self.prev_dq[0]-self.prev_dq[1])


    def solve(self):
        self.calc.start_timing('solve')
        # TODO: enable fixed dq-calculations in SCC (for band-structures)
        # TODO: move this from here to elements
        #self.calc.ia.check_too_close_distances() 
        dq=self.guess_dq()
        self.H0, self.S, self.dH0, self.dS = self.calc.ia.get_matrices()
#        try:
        if True:
            if self.SCC:
                raise NotImplementedError
                self.es.construct_Gamma_matrix()
            self.e, self.wf = self.solver.get_states(self.calc,dq,self.H0,self.S,self.count)
            self.calc.el.set_solved('ground state')
            #self.check_mulliken_charges()
            self.large_update()
            self.count+=1
            self.calc.stop_timing('solve')
#        except Exception, ex:
#            self.calc.stop_timing('solve')
#            raise Exception(ex)


    def check_mulliken_charges(self):
        """ Check that the Mulliken populations are physically
        reasonable. """
        dQ = self.mulliken()
        if self.calc.verbose_SCC:
            print "Mulliken populations: min=%0.3f, max=%0.3f" % (nu.min(dQ),nu.max(dQ))
        Z = self.calc.el.get_atomic_numbers()
        for dq, z in zip(dQ, Z):
            if dq < -z or dq > z:
                for dq, z in zip(dQ, Z):
                    print >> self.calc.get_output(), "Z=%i    dq=%0.3f    excess charge=%0.3f" % (z, dq, -dq)
                raise Exception("The Mulliken charges are insane!")


    def update(self,e,wf):
        """ Update all essential stuff from given wave functions. """
        self.calc.start_timing('update')
        self.e=e
        self.wf=wf
        self.f=self.occu.occupy(e)
        self.calc.start_timing('rho')
        # TODO: rhoS etc...
#        self.rho0=fortran_rho0(self.wf,self.f,self.calc.el.nr_ia_orbitals,self.calc.el.ia_orbitals,self.norb)
        self.rho = fortran_rhoc(self.wf,self.f,self.norb,self.nk)
        
        
#        for ik in range(self.nk):
#            norm=0.0
#            a=3
#            print nu.dot( nu.dot(self.wf[ik,:,:].transpose().conjugate(),self.S[ik,:,:]),self.wf[ik,:,:] )
#            for p in range(self.norb):
#                norm+=nu.dot(self.wf[ik,:,a].conjugate(),self.S[ik,:,p])*self.wf[ik,a,p]
#            print 'norm',norm    
#            rr = self.rho[ik,:,:]
#            ss = self.S[ik,:,:]
#            print 'trace',ik,2*(nu.trace(nu.dot(rr,ss)) + nu.trace(nu.dot(ss,rr))).real 
        self.calc.stop_timing('rho')
#        self.rho0S_diagonal=matmul_diagonal(self.rho0,self.S,self.norb)
        if self.SCC:
            raise NotImplementedError
            self.dq=self.mulliken()
            self.es.set_dq(self.dq)
        self.calc.stop_timing('update')


    def large_update(self):
        """ Update stuff from eigenstates needed later for forces etc. """
        self.calc.start_timing('final update')
        self.Swf = None
        #self.rho=fortran_rho(self.wf,self.f,self.norb) # the complete density matrix (even needed?)
        self.dH = nu.zeros_like(self.dH0)
        if self.SCC:
            raise NotImplementedError
            self.prev_dq=[self.dq, self.prev_dq[0]]
            for a in range(3):
                self.dH[:,:,a]=self.dH0[:,:,a] + self.es.get_h1()*self.dS[:,:,a]
        else:
            self.dH = self.dH0

        # density matrix weighted by eigenenergies
        # TODO: rhoe0 ...
#        self.rhoe0=fortran_rhoe0(self.wf,self.f,self.e,self.calc.el.nr_ia_orbitals,self.calc.el.ia_orbitals,self.norb)
        self.rhoe=fortran_rhoec(self.wf,self.f,self.e,self.norb,self.nk)
        self.calc.stop_timing('final update')

    def get_dq(self):
        return self.dq

    def get_eigenvalues(self):
        return self.e

    def get_occupations(self):
        return self.f

    def get_homo(self):
        """ Return highest orbital with occu>0.99. """
        raise NotImplementedError
        for i in range(self.norb)[::-1]:
            if self.f[i]>0.99: return i

    def get_lumo(self):
        """ Return lowest orbital with occu<1.01. """
        raise NotImplementedError
        for i in range(self.norb):
            if self.f[i]<1.01: return i

    def get_hoc(self):
        """ Return highest partially occupied orbital. """
        raise NotImplementedError
        for i in range(self.norb)[::-1]:
            if self.f[i]>1E-9: return i


    def get_fermi_level(self):
        raise NotImplementedError('getting Fermi-level is a matter of occu-instance. Use it instead.')


    def mulliken(self):
        """ Return excess Mulliken populations. """
        raise NotImplementedError
        q=[]
        for i in range(self.nat):
            orbitals=self.calc.el.orbitals(i,indices=True)
            q.append( sum(self.rho0S_diagonal[orbitals]) )
        return nu.array(q)-self.calc.el.get_valences()


    def mulliken_transfer(self,k,l):
        """ Return Mulliken transfer charges between states k and l. """
        raise NotImplementedError
        if self.Swf==None:
            self.Swf=symmetric_matmul(self.S,self.wf)
        q=[]
        for i in range(self.nat):
            iorb=self.calc.el.orbitals(i,indices=True)
            qi=sum( [self.wf[oi,k]*self.Swf[oi,l]+self.wf[oi,l]*self.Swf[oi,k] for oi in iorb] )
            q.append(qi/2)
        return nu.array(q)


    def band_structure_energy(self):
        """ Return band structure energy. """
        # TODO: use rho0-version
        ebs=0.0
        for ik in range(self.nk):
            ebs += self.wk[ik] * nu.trace(nu.dot(self.rho[ik],self.H0[ik]))
        assert ebs.imag<1E-13
        
#        eb2=0.0
#        for ik in range(self.nk):
#            for a in range(self.norb):
#                eb2+=self.e[ik,a]*self.wk[ik]*self.f[ik,a]
#        print 'ebs,eb2',ebs,eb2
#        
#        
        return ebs.real 


    def band_structure_forces(self):
        """ Return forces from band structure. """
        self.calc.start_timing('f_bs')
        norbs=self.calc.el.nr_orbitals
        inds=self.calc.el.atom_orb_indices2
        
        f=fortran_fbsc(self.rho,self.rhoe,self.dH,self.dS,norbs,inds,self.wk,self.norb,self.nat,self.nk)
        self.calc.stop_timing('f_bs')
        return f


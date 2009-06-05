# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

from solver import Solver
from electrostatics import Electrostatics
from occupations import Occupations
from hotbit.fortran.misc import fortran_rho
from hotbit.fortran.misc import fortran_rho0
from hotbit.fortran.misc import fortran_rhoe0
from hotbit.fortran.misc import symmetric_matmul
from hotbit.fortran.misc import matmul_diagonal
from hotbit.fortran.misc import fortran_fbs
import numpy as nu
from weakref import proxy

class States:

    def __init__(self,calc):
        self.es=Electrostatics(calc)
        self.solver=Solver(calc)
        self.calc=proxy(calc)
        width=calc.get('width')
        self.occu=Occupations(calc.el.get_number_of_electrons(),width=width)
        self.nat=len(calc.el)
        self.norb=calc.el.get_nr_orbitals()
        self.prev_dq=[None,None]
        self.count=0
        self.SCC=calc.get('SCC')
        self.rho=None
        self.rhoe0=None

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
        self.calc.ia.check_too_close_distances() #TODO: move this from here to elements
        dq=self.guess_dq()
        self.H0, self.S, self.dH0, self.dS = self.calc.ia.get_matrices()
        try:
            if self.SCC:
                self.es.construct_Gamma_matrix()
            self.e, self.wf = self.solver.get_states(self.calc,dq,self.H0,self.S,self.count)
            self.calc.el.set_solved('ground state')
            self.check_mulliken_charges()
            self.large_update()
            self.count+=1
            self.calc.stop_timing('solve')
        except Exception, ex:
            self.calc.stop_timing('solve')
            raise Exception(ex)


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
        self.rho0=fortran_rho0(self.wf,self.f,self.calc.el.nr_ia_orbitals,self.calc.el.ia_orbitals,self.norb)
        self.calc.stop_timing('rho')
        self.rho0S_diagonal=matmul_diagonal(self.rho0,self.S,self.norb)
        if self.SCC:
            self.dq=self.mulliken()
            self.es.set_dq(self.dq)
        self.calc.stop_timing('update')


    def large_update(self):
        """ Update stuff from eigenstates needed later for forces etc. """
        self.calc.start_timing('final update')
        self.Swf=None
        #self.rho=fortran_rho(self.wf,self.f,self.norb) # the complete density matrix (even needed?)
        self.dH=nu.zeros_like(self.dH0)
        if self.SCC:
            self.prev_dq=[self.dq, self.prev_dq[0]]
            for a in range(3):
                self.dH[:,:,a]=self.dH0[:,:,a] + self.es.get_h1()*self.dS[:,:,a]
        else:
            self.dH=self.dH0

        # density matrix weighted by eigenenergies
        self.rhoe0=fortran_rhoe0(self.wf,self.f,self.e,self.calc.el.nr_ia_orbitals,self.calc.el.ia_orbitals,self.norb)
        self.calc.stop_timing('final update')

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

    def get_fermi_level(self):
        """ Return the Fermi level of the system. The Fermi level
            is defined by the fit of 

                               2
              f(e) =  --------------------
                      exp((e-e_F)/kbT) + 1

            to the eigenenergy-occupation -points (kbT is the width of
            the Fermi distribution and e_F the Fermi energy).
            If the kbT is smaller than 0.004 Ha the eigenenergy
            of the HOMO is considered as the Fermi level. """
        e = self.get_eigenvalues()
        f = self.get_occupations()
        from scipy.optimize import leastsq
        fermifunc = lambda p, e:    2./(nu.exp((e-p[0])/p[1]) + 1)
        errfunc   = lambda p, e, y: fermifunc(p, e) - y
        # p = [Fermi_level, width]
        p = [0.5 * (e[self.get_homo()] + e[self.get_lumo()]), 1.0]
        p, success = leastsq(errfunc, p[:], args=(e,f), ftol=0.0001, xtol=0.001)
        e_g = nu.linspace(e[0], e[-1], 500)
        #import pylab
        #pylab.figure()
        #pylab.plot(e,f,"b^", e_g, fermifunc(p, e_g), "r-")
        #pylab.show()
        if success:
            if p[1] < 0.004:
                return e[self.get_homo()]
            else:
                return p[0]
        else:
            raise RuntimeError("Could not define the Fermi level.")


    def mulliken(self):
        """ Return excess Mulliken populations. """
        q=[]
        for i in range(self.nat):
            orbitals=self.calc.el.orbitals(i,indices=True)
            q.append( sum(self.rho0S_diagonal[orbitals]) )
        return nu.array(q)-self.calc.el.get_valences()

    def mulliken_transfer(self,k,l):
        """ Return Mulliken transfer charges between states k and l. """
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
        return nu.trace(nu.dot(self.rho0,self.H0))

    def band_structure_forces(self):
        """ Return forces arising from band structure. """
        self.calc.start_timing('f_bs')

        norbs=self.calc.el.nr_orbitals
        inds=self.calc.el.atom_orb_indices2

        f=fortran_fbs(self.rho0,self.rhoe0,self.dH,self.dS,norbs,inds,self.norb,self.nat)
        self.calc.stop_timing('f_bs')
        return -2*nu.real(f)


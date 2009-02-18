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
            return nu.zeros((n,)) - self.calc.get('charge')/n
        elif self.count==1:    # use previous charges
            return self.prev_dq[0]
        else:                  # linear extrapolation
            return self.prev_dq[0] + (self.prev_dq[0]-self.prev_dq[1])


    def solve(self):
        self.calc.start_timing('solve')
        self.calc.ia.check_too_close_distances()
        dq=self.guess_dq()
        self.H0, self.S, self.dH0, self.dS=self.calc.ia.get_matrices()
        try:
            if self.SCC:
                self.es.construct_tables()
            self.e, self.wf=self.solver.get_states(self.calc,dq,self.H0,self.S,self.count)
            self.calc.el.set_solved('ground state')
            self.large_update()
            self.count+=1
            self.calc.stop_timing('solve')
        except Exception, ex:
            self.calc.stop_timing('solve')
            raise Exception(ex)


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
                self.dH[:,:,a]=self.dH0[:,:,a] + self.es.get_H1()*self.dS[:,:,a]
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

    def mulliken(self):
        """ Return excess Mulliken populations. """
        q=[]
        for i in range(self.nat):
            orbitals=self.calc.el.orbitals(i,indices=True)
            q.append( sum(self.rho0S_diagonal[orbitals]) )
        return nu.array(q)-self.calc.el.get_valences()

    def mulliken_I(self, I):
        """ Return the Mulliken population on atom I. """
        rho_tilde = 0.5*(self.rho0 + self.rho0.conjugate().transpose())
        rho_tilde_S = nu.dot(rho_tilde, self.S)
        q = self.trace_I(I, rho_tilde_S)
        return q

    def mulliken_mu(self, mu):
        """ Return the population of basis state mu. """
        raise NotImplementedError('Check what this means')
        rho_tilde = 0.5*(self.rho0 + self.rho0.conjugate().transpose())
        rho_tilde_S = nu.dot(rho_tilde, self.S)
        return rho_tilde_S[mu,mu]

    def mulliken_I_k(self, I, k):
        """ Return the Mulliken population of atom I from eigenstate k. """
        rho_k = self.get_rho_k(k)
        rho_tilde_k = 0.5*(rho_k + rho_k.conjugate().transpose())
        q_Ik = self.trace_I(I, nu.dot(rho_tilde_k,self.S))
        return q_Ik

    def mulliken_I_l(self, I, l):
        """ Return the Mulliken population of atom I basis states with
        angular momentum l. """
        rho_tilde = 0.5*(self.rho0 + self.rho0.conjugate().transpose())
        rho_tilde_S = nu.dot(rho_tilde, self.S)
        orb_indices = self.calc.el.orbitals(I, indices=True)
        orbs = self.calc.el.orbitals(I)
        q = 0
        for i, orb in zip(orb_indices, orbs):
            if   's' in orb['orbital']: l_orb = 0
            elif 'p' in orb['orbital']: l_orb = 1
            elif 'd' in orb['orbital']: l_orb = 2
            else: raise RuntimeError('Something wrong with orbital types')
            if l_orb == l:
                q += rho_tilde_S[i,i]
        return q

    def mulliken_I_k_l(self, I, k, l):
        """ Return the Mulliken population of atom I from eigenstate k
            from basis functions with angular momentum l. """
        rho_k = self.get_rho_k(k)
        rho_tilde_k = 0.5*(rho_k + rho_k.conjugate().transpose())
        rho_tilde_k_S = nu.dot(rho_k, self.S)
        q_Ikl = 0.0
        orb_indices = self.calc.el.orbitals(I, indices=True)
        orbs = self.calc.el.orbitals(I)
        for i, orb in zip(orb_indices, orbs):
            if   's' in orb['orbital']: l_orb = 0
            elif 'p' in orb['orbital']: l_orb = 1
            elif 'd' in orb['orbital']: l_orb = 2
            else: raise RuntimeError('Something wrong with orbital types')
            if l_orb == l:
                q_Ikl += rho_tilde_k_S[i,i]
        return q_Ikl

    def get_rho_k(self, k):
        rho_k = nu.zeros_like(self.rho0)
        for i, c_ik in enumerate(self.wf[:,k]):
            for j, c_jk in enumerate(self.wf[:,k]):
                rho_k[i,j] = c_ik*c_jk.conjugate()
        return rho_k

    def trace_I(self, I, matrix):
        """ Return partial trace over atom I's orbitals. """
        ret = 0
        I = self.calc.el.orbitals(I, indices=True)
        for i in I:
            ret += matrix[i,i]
        return ret

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


    def mayer_bond_order(self, a, b):
        """ Returns the Mayer bond-order of the bond between the
        atoms A and B (a and b are atom indices). """
        assert type(a) == int
        assert type(b) == int
        A = self.calc.el.orbitals(a, indices=True)
        B = self.calc.el.orbitals(b, indices=True)
        rho_tilde = 0.5*(self.rho0 + self.rho0.conjugate().transpose())
        rho_tilde_S = nu.dot(rho_tilde, self.S)
        B_AB = 0
        for i in A:
            for j in B:
                B_AB += rho_tilde_S[i,j]*rho_tilde_S[j,i]
        return B_AB

    def get_covalent_energy(self):
        """ Returns the covalent bond energy of the whole system. """
        E_bs = self.band_structure_energy()
        rho = self.rho0
        H0 = self.H0
        S = self.S
        epsilon = nu.zeros_like(H0)
        for i in range(len(H0)):
            for j in range(len(H0)):
                epsilon[i,j] = 0.5*(H0[i,i] + H0[j,j])
        E = nu.sum(rho*epsilon*S)
        return E_bs - E

    def E_cov_m_n(self, m, n, sigma, npts=400):
        """ Return the covalent energy of the orbital pairs mu and nu
        as a function of an energy. Returns the energy grid and covalent
        energy grid. """
        eigs = self.get_eigenvalues()
        e_min, e_max = min(eigs), max(eigs)
        e_range = e_max - e_min
        e_min -= e_range*0.1
        e_max += e_range*0.1
        e_g = nu.linspace(e_min, e_max, npts)
        f = self.get_occupations()
        E_cov_mn = nu.zeros_like(e_g)
        epsilon = nu.zeros_like(self.H0)
        for i in range(len(self.H0)):
            for j in range(len(self.H0)):
                epsilon[i,j] = 0.5*(self.H0[i,i] + self.H0[j,j])
        for k, e_k in enumerate(eigs):
            f_k = f[k]
            rho_k = self.get_rho_k(k)
            for i, e in enumerate(e_g):
                E_cov_mn[i] += f_k*nu.exp(-(e-e_k)**2/(2*sigma**2))*rho_k[m,n]*(self.H0[n,m]+epsilon[n,m]*self.S[n,m])
        return e_g, E_cov_mn

    def E_cov_mn(self):
        """ Return the covalent energy of the orbital pairs mu and nu
        as a function of an energy. """
        f = self.get_occupations()
        E_cov_mn = 0.0
        epsilon = nu.zeros_like(self.H0)
        for i in range(len(self.H0)):
            for j in range(len(self.H0)):
                epsilon[i,j] = 0.5*(self.H0[i,i] + self.H0[j,j])
        for k, f_k in enumerate(f):
            rho_k = self.get_rho_k(k)
            E_cov_mn += f_k*nu.trace(nu.dot(rho_k, self.H0)) - f_k*nu.sum(rho_k*epsilon*self.S)
        return E_cov_mn

    def DOS(self, sigma, npts=300):
        """ Return the energy array and corresponding DOS array. """
        eigs = self.get_eigenvalues()
        e_min, e_max = min(eigs), max(eigs)
        e_range = e_max - e_min
        e_min -= e_range*0.1
        e_max += e_range*0.1
        e_g = nu.linspace(e_min, e_max, npts)
        ldos = nu.zeros_like(e_g)
        for e_k in eigs:
            ldos += [nu.exp(-(e-e_k)**2/(2*sigma**2)) for e in e_g]
        return e_g, ldos

    def LDOS(self, sigma, indices=None, npts=300):
        """ Return the energy array and corresponding LDOS array
            calculated using Mulliken population analysis. Indices
            refer to the atoms that are included to the LDOS. """
        eigs = self.get_eigenvalues()
        e_min, e_max = min(eigs), max(eigs)
        e_range = e_max - e_min
        e_min -= e_range*0.1
        e_max += e_range*0.1
        e_g = nu.linspace(e_min, e_max, npts)
        ldos = nu.zeros_like(e_g)
        if indices == None:
            indices = range(self.nat)
        elif type(indices) == int:
            indices = [indices]
        N_el = self.calc.el.get_number_of_electrons()
        for I in indices:
            for k, e_k in enumerate(eigs):
                q_Ik = self.mulliken_I_k(I,k)
                ldos_k = [nu.exp(-(e-e_k)**2/(2*sigma**2)) for e in e_g]
                ldos += nu.array(ldos_k) * q_Ik
        return e_g, ldos

    def PDOS(self, sigma, indices=None, l='spd', npts=300):
        """ Return the energy array and corresponding PDOS array
            calculated using Mulliken population analysis. Indices refer
            to the atoms and l to the angular momenta that are included. """
        eigs = self.get_eigenvalues()
        e_min, e_max = min(eigs), max(eigs)
        e_range = e_max - e_min
        e_min -= e_range*0.1
        e_max += e_range*0.1
        e_g = nu.linspace(e_min, e_max, npts)
        pdos = nu.zeros_like(e_g)
        if indices == None:
            indices = range(self.nat)
        elif type(indices) == int:
            indices = [indices]
        if l == 'spd':
            l = [0,1,2]
        elif type(l) == str:
            l = self.get_angular_momenta(l)
        else:
            raise RuntimeError("l must be orbital types, for example l='sp'")
        for li in l:
            for k, e_k in enumerate(eigs):
                pdos_k = [nu.exp(-(e-e_k)**2/(2*sigma**2)) for e in e_g]
                for I in indices:
                    q_Ikl = self.mulliken_I_k_l(I,k,li)
                    pdos += nu.array(pdos_k) * q_Ikl
        return e_g, pdos

    def get_angular_momenta(self, l):
        ret = []
        if 's' in l: ret.append(0)
        if 'p' in l: ret.append(1)
        if 'd' in l: ret.append(2)
        return ret

    def hybridization(self, la, lb):
        """ Return a number between zero and one that describes how much the
            wave functions of the atoms are hybridized between angular
            momenta la and lb. """
        h = 0.0
        for k, fk in enumerate(self.f):
            for I in range(len(self.calc.el)):
                h += fk * self.mulliken_I_k_l(I, k, la) * self.mulliken_I_k_l(I, k, lb)
        return h


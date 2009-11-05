from ase import Hartree
import numpy as nu


def get_angular_momenta(l):
    ret = []
    if 's' in l: ret.append(0)
    if 'p' in l: ret.append(1)
    if 'd' in l: ret.append(2)
    return ret


class MullikenAnalysis:
    """ A class that calculates different Mulliken charges. All the units
    are in electronvolts and Angstroms. """

    def __init__(self, calc):
        self.calc = calc
        self.rho = calc.st.rho
        self.H0 = calc.st.H0*Hartree
        self.S = calc.st.S

        self.rho_tilde = 0.5*(self.rho + self.rho.conjugate().transpose())
        self.rho_tilde_S = nu.dot(self.rho_tilde, self.S)
        self.eigs = self.calc.st.get_eigenvalues().copy()*Hartree
        if calc.st.get_lumo() == None:
            fermi_energy = self.eigs[calc.st.get_homo()]
        else:
            fermi_energy = 0.5*(self.eigs[calc.st.get_homo()] + self.eigs[calc.st.get_lumo()])
        self.eigs -= fermi_energy


    def delta_sigma(self, x, x0, sigma):
        """ Return the value of a gaussian function centered at x0
        with variance sigma^2 at point x. """
        return 1./nu.sqrt(2*nu.pi*sigma**2)*nu.exp(-(x-x0)**2/(2*sigma**2))


    def get_rho_k(self, k):
        """ Return the k-contribution of the density matrix without
        multiplication with the occupation number f_k. """
        return nu.outer(self.calc.st.wf[:,k],self.calc.st.wf[:,k].conjugate())


    def trace_I(self, I, matrix):
        """ Return partial trace over atom I's orbitals. """
        ret = 0
        I = self.calc.el.orbitals(I, indices=True)
        for i in I:
            ret += matrix[i,i]
        return ret


    def mulliken(self):
        """ Return excess Mulliken populations. """
        q=[]
        for i in range(len(self.calc.el)):
            q.append( nu.sum(self.trace_I(i, self.rho_tilde_S)) )
        return nu.array(q)-self.calc.el.get_valences()


    def mulliken_I(self, I):
        """ Return the Mulliken population on atom I. """
        return self.trace_I(I, self.rho_tilde_S)


    def mulliken_mu(self, mu):
        """ Return the population of basis state mu. """
        return self.rho_tilde_S[mu,mu]


    def mulliken_I_k(self, I, k):
        """ Return the Mulliken population of atom I from eigenstate k. """
        rho_k = self.get_rho_k(k)
        rho_tilde_k = 0.5*(rho_k + rho_k.conjugate().transpose())
        q_Ik = self.trace_I(I, nu.dot(rho_tilde_k,self.S))
        return q_Ik


    def mulliken_I_l(self, I, l):
        """ Return the Mulliken population of atom I basis states with
        angular momentum l. """
        orb_indices = self.calc.el.orbitals(I, indices=True)
        orbs = self.calc.el.orbitals(I)
        q = 0
        for i, orb in zip(orb_indices, orbs):
            if   's' in orb['orbital']: l_orb = 0
            elif 'p' in orb['orbital']: l_orb = 1
            elif 'd' in orb['orbital']: l_orb = 2
            else: raise RuntimeError('Something wrong with orbital types')
            if l_orb == l:
                q += self.rho_tilde_S[i,i]
        return q


    def mulliken_I_k_l(self, I, k, l):
        """ Return the Mulliken population of atom I from eigenstate k
            from basis functions with angular momentum l. """
        rho_k = self.get_rho_k(k)
        rho_tilde_k = 0.5*(rho_k + rho_k.conjugate().transpose())
        rho_tilde_k_S = nu.dot(rho_tilde_k, self.S)
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


class MullikenBondAnalysis(MullikenAnalysis):
    """ A class that calculates analyzes atom bonds using the Mulliken
    charges. """


    def __init__(self, calc):
        MullikenAnalysis.__init__(self, calc)
        self.bar_epsilon = nu.zeros_like(self.H0)
        for i in range(len(self.H0)):
            for j in range(len(self.H0)):
                self.bar_epsilon[i,j] = 0.5*(self.H0[i,i] + self.H0[j,j])


    def get_mayer_bond_order(self, a, b):
        """ Returns the Mayer bond-order of the bond between the
        atoms A and B (a and b are atom indices). """
        assert type(a) == int
        assert type(b) == int
        A = self.calc.el.orbitals(a, indices=True)
        B = self.calc.el.orbitals(b, indices=True)
        B_AB = 0
        for i in A:
            for j in B:
                B_AB += self.rho_tilde_S[i,j]*self.rho_tilde_S[j,i]
        return B_AB


    def get_covalent_energy(self):
        """ Returns the covalent bond energy of the whole system. """
        E_bs = self.calc.st.band_structure_energy()*Hartree
        E = nu.sum(self.rho*self.bar_epsilon*self.S)
        return E_bs - E


    def get_E_cov_m_n(self, m, n, sigma, npts=500, e_min=None, e_max=None, occupations=False):
        """ Return the covalent energy of the orbital pairs mu and nu
        as a function of an energy. Returns the energy grid and covalent
        energy grid, where the fermi-energy is shifted to zero. """
        if e_min == None:
            e_min = min(self.eigs)
        if e_max == None:
            e_max = max(self.eigs)
            e_range = e_max - e_min
            e_min -= e_range*0.1
            e_max += e_range*0.1
        e_g = nu.linspace(e_min, e_max, npts)
        f = self.calc.st.get_occupations()
        E_cov_mn = nu.zeros_like(e_g)
        mat_nm = self.H0[n,m] - self.bar_epsilon[n,m]*self.S[n,m]
        for k, e_k in enumerate(self.eigs):
            rho_k_mn = self.calc.st.wf[m,k]*self.calc.st.wf[n,k].conjugate()
            if occupations:
                rho_k_mn *= f[k]
            E_cov_k = rho_k_mn*mat_nm
            add = self.delta_sigma(e_g, e_k, sigma)*E_cov_k
            E_cov_mn += add
        return e_g, E_cov_mn


    def get_E_cov_I_J(self, I, J, sigma, npts=500, e_min=None, e_max=None, occupations=False):
        """ Returns the energy grid and corresponding covalent bonding
            energy values for the atom pair I-J. """
        I = self.calc.el.orbitals(I, indices=True)
        J = self.calc.el.orbitals(J, indices=True)
        E_cov_IJ = nu.zeros(npts)
        e_g = None
        for m in I:
            for n in J:
                e_g, E_cov_mn = self.get_E_cov_m_n(m, n, sigma, npts, e_min, e_max, occupations=occupations)
                E_cov_IJ += E_cov_mn
        return e_g, 2*E_cov_IJ


    def get_E_cov_la_lb(self, la, lb, sigma, npts=500, e_min=None, e_max=None, occupations=False):
        print "Covalent energy graph for angular momenta %s and %s..." % (la, lb)
        e_g = None
        E_cov_lalb = nu.zeros(npts)
        if len(la) != 1 or len(lb) != 1:
            raise Exception('Give only one angular momentum at a time.')
        if e_min == None:
            e_min = min(self.eigs)
        if e_max == None:
            e_max = max(self.eigs)
            e_range = e_max - e_min
            e_min -= e_range*0.1
            e_max += e_range*0.1
        e_g = nu.linspace(e_min, e_max, npts)
        E_cov_mn = nu.zeros_like(e_g)
        f = self.calc.st.get_occupations()
        for k, e_k in enumerate(self.eigs):
            print "  %s%s: Analyzing state %i/%i" % (la, lb, k+1, len(self.eigs))
            rho_k = self.get_rho_k(k)
            if occupations:
                rho_k *= f[k]
            mat = self.H0 - self.bar_epsilon*self.S
            add = self.delta_sigma(e_g, e_k, sigma)
            w_k = 0.0
            for i, orb_i in enumerate(self.calc.el.orb):
                for j, orb_j in enumerate(self.calc.el.orb):
                    orb_type_i = orb_i['orbital']
                    orb_type_j = orb_j['orbital']
                    if (la == orb_type_i[0] and lb == orb_type_j[0]) or\
                       (la == orb_type_j[0] and lb == orb_type_i[0]):
                        w_k += rho_k[i,j]*mat[j,i]
            E_cov_mn += add*w_k
        return e_g, E_cov_mn


    def A_I(self, I):
        """ Return the absolute energy of atom I. """
        gamma_II = self.calc.st.es.gamma(I,I)*Hartree
        dq_I = self.mulliken_I(I) - self.calc.el.get_valences()[I]
        return 0.5*gamma_II*dq_I**2 + self.E_prom_I(I)


    def E_prom_I(self, I):
        """ Return the promotion energy of atom I. """
        orb_indices = self.calc.el.orbitals(I, indices=True)
        ret = 0.0
        for m in orb_indices:
            q_mu = self.mulliken_mu(m)
            q_mu_free = self.calc.el.get_free_population(m)
            ret += (q_mu-q_mu_free)*self.H0[m,m]
        return ret


    def B_IJ(self, I, J):
        """ Return the absolute bond energy between atoms I and J. """
        if I == J:
             return 0
        dist_IJ = self.calc.el.distance(I, J)

        sI = self.calc.el.symbol(I)
        sJ = self.calc.el.symbol(J)
        V_rep_IJ = self.calc.rep.vrep[sI+sJ](dist_IJ)*Hartree

        gamma_IJ = self.calc.st.es.gamma(I, J)*Hartree
        dq_I = self.mulliken_I(I) - self.calc.el.get_valences()[I]
        dq_J = self.mulliken_I(J) - self.calc.el.get_valences()[J]

        ret = V_rep_IJ + gamma_IJ*dq_I*dq_J
        for m in self.calc.el.orbitals(I, indices=True):
            for n in self.calc.el.orbitals(J, indices=True):
                mat = self.rho[n,m]*self.H0[m,n]-self.rho_tilde[n,m]*self.S[m,n]*self.bar_epsilon[m,n]
                ret += (mat + mat.conjugate())
        return ret


    def get_AB_I(self, I):
        """ Return the atom I's contribution to the total binding
        energy of the system. """
        ret = self.A_I(I)
        for J in range(len(self.calc.el)):
            if J != I:
                ret += 0.5*self.B_IJ(I,J)
        return ret


class DensityOfStates(MullikenAnalysis):
    """ A class that calculates different kinds of local and projected
    density of states using the Mulliken charges. """

    def __init__(self, calc):
        MullikenAnalysis.__init__(self, calc)


    def DOS(self, sigma, npts=500, e_min=None, e_max=None, occupations=False):
        """ Return the energy array and corresponding DOS array.
            
            sigma: the broadening of the Gaussian distribution
            npts:  the number of grid points
            e_min: the minimum of the energy array
            e_max: the maximum of the energy array
            occupations: weight the states with the occupation number of
                         the state.
        """
        if e_min == None:
            e_min = min(self.eigs)
        if e_max == None:
            e_max = max(self.eigs)
        e_range = e_max - e_min
        e_min -= e_range*0.1
        e_max += e_range*0.1
        e_g = nu.linspace(e_min, e_max, npts)
        dos = nu.zeros_like(e_g)
        f = self.calc.st.f
        for k, e_k in enumerate(self.eigs):
            add = nu.array([self.delta_sigma(e, e_k, sigma) for e in e_g])
            if occupations:
                add *= f[k]
            dos += add
        return e_g, dos


    def LDOS(self, sigma, indices=None, npts=500, e_min=None, e_max=None, occupations=False):
        """ Return the energy array and corresponding local DOS array
            calculated using Mulliken population analysis.
            
            sigma: the broadening of the Gaussian distribution
            indices: the indices of atoms included to the LDOS
            npts:  the number of grid points
            e_min: the minimum of the energy array
            e_max: the maximum of the energy array
            occupations: weight the states with the occupation number of
                         the state.
        """
        if e_min == None:
            e_min = min(self.eigs)
        if e_max == None:
            e_max = max(self.eigs)
        e_range = e_max - e_min
        e_min -= e_range*0.1
        e_max += e_range*0.1
        e_g = nu.linspace(e_min, e_max, npts)
        ldos = nu.zeros_like(e_g)
        if indices == None:
            indices = range(len(self.calc.el))
        elif type(indices) == int:
            indices = [indices]
        N_el = self.calc.el.get_number_of_electrons()
        f = self.calc.st.f
        for I in indices:
            for k, e_k in enumerate(self.eigs):
                q_Ik = self.mulliken_I_k(I,k)
                ldos_k = [self.delta_sigma(e, e_k, sigma) for e in e_g]
                if occupations:
                    ldos += nu.array(ldos_k) * q_Ik * f[k]
                else:
                    ldos += nu.array(ldos_k) * q_Ik
        return e_g, ldos


    def PDOS(self, sigma, indices=None, l='spd', npts=500, e_min=None, e_max=None, occupations=False):
        """ Return the energy array and corresponding projected DOS array
            calculated using Mulliken population analysis. Indices refer
            to the atoms and l to the angular momenta that are included.
            
            sigma: the broadening of the Gaussian distribution
            indices: the indices of atoms included to the LDOS
            npts:  the number of grid points
            e_min: the minimum of the energy array
            e_max: the maximum of the energy array
            occupations: weight the states with the occupation number of
                         the state.
        """
        if e_min == None:
            e_min = min(self.eigs)
        if e_max == None:
            e_max = max(self.eigs)
        e_range = e_max - e_min
        e_min -= e_range*0.1
        e_max += e_range*0.1
        e_g = nu.linspace(e_min, e_max, npts)
        pdos = nu.zeros_like(e_g)
        if indices == None:
            indices = range(len(self.calc.el))
        elif type(indices) == int:
            indices = [indices]
        if l == 'spd':
            l = [0,1,2]
        elif type(l) == str:
            l = get_angular_momenta(l)
        else:
            raise RuntimeError("l must be orbital types, for example l='sp'")
        for li in l:
            f = self.calc.st.f
            for k, e_k in enumerate(self.eigs):
                pdos_k = nu.array([self.delta_sigma(e, e_k, sigma) for e in e_g])
                if occupations:
                    pdos_k *= f[k]
                for I in indices:
                    q_Ikl = self.mulliken_I_k_l(I,k,li)
                    pdos += pdos_k * q_Ikl
        return e_g, pdos




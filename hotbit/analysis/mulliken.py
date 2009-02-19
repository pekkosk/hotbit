import numpy as nu


def get_angular_momenta(l):
    ret = []
    if 's' in l: ret.append(0)
    if 'p' in l: ret.append(1)
    if 'd' in l: ret.append(2)
    return ret


class MullikenAnalysis:
    """ A class that calculates different Mulliken charges. """

    def __init__(self, calc):
        self.calc = calc
        self.rho = calc.st.rho0
        self.H0 = calc.st.H0
        self.S = calc.st.S

        self.rho_tilde = 0.5*(self.rho + self.rho.conjugate().transpose())
        self.rho_tilde_S = nu.dot(self.rho_tilde, self.S)
        self.rho_k = nu.zeros((len(self.calc.st.get_eigenvalues()), len(self.rho), len(self.rho)))
        self.built_rho_k = [False for i in range(len(self.rho_k))]



    def get_rho_k(self, k):
        """ Return the k-contribution of the density matrix without
        multiplication with the occupation number f_k. """
        if self.built_rho_k[k]:
            return self.rho_k[k]
        else:
            for i, c_ik in enumerate(self.calc.st.wf[:,k]):
                for j, c_jk in enumerate(self.calc.st.wf[:,k]):
                    self.rho_k[k,i,j] = c_ik*c_jk.conjugate()
            self.built_rho_k[k] = True
        return self.rho_k[k]


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
        raise NotImplementedError('Check what this means')
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


class MullikenBondAnalysis(MullikenAnalysis):
    """ A class that calculates analyzes atom bonds using the Mulliken
    charges. """

    def __init__(self, calc):
        MullikenAnalysis.__init__(self, calc)
        self.calc = calc
        self.bar_epsilon = nu.zeros_like(self.H0)
        for i in range(len(self.H0)):
            for j in range(len(self.H0)):
                self.bar_epsilon[i,j] = 0.5*(self.H0[i,i] + self.H0[j,j])


    def hybridization(self, la, lb):
        """ Return a number between zero and one that describes how much the
            wave functions of the atoms are hybridized between angular
            momenta la and lb. """
        h = 0.0
        for k, fk in enumerate(self.calc.st.get_occupations()):
            for I in range(len(self.calc.el)):
                h += fk * self.mulliken_I_k_l(I, k, la) * self.mulliken_I_k_l(I, k, lb)
        return h


    def delta_sigma(self, x, x0, sigma):
        """ Return the value of a gaussian function centered at x0
        with variance sigma^2 at point x. """
        return nu.sqrt(2*nu.pi*sigma**2)*nu.exp(-(x-x0)**2/(2*sigma**2))


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
        E_bs = self.calc.st.band_structure_energy()
        E = nu.sum(self.rho*self.bar_epsilon*self.S)
        return E_bs - E


    def E_cov_m_n(self, m, n, sigma, npts=800):
        """ Return the covalent energy of the orbital pairs mu and nu
        as a function of an energy. Returns the energy grid and covalent
        energy grid. """
        eigs = self.calc.st.get_eigenvalues()
        e_min, e_max = min(eigs), max(eigs)
        e_range = e_max - e_min
        e_min -= e_range*0.1
        e_max += e_range*0.1
        e_g = nu.linspace(e_min, e_max, npts)
        f = self.calc.st.get_occupations()
        E_cov_mn = nu.zeros_like(e_g)
        mat = self.H0[n,m] - self.bar_epsilon[n,m]*self.S[n,m]
        for k, e_k in enumerate(eigs):
            f_k = f[k]
            rho_k = self.get_rho_k(k)
            E_cov_k = rho_k[m,n]*mat
            for i, e in enumerate(e_g):
                E_cov_mn[i] += f_k/nu.sqrt(2*nu.pi*sigma**2)*nu.exp(-(e-e_k)**2/(2*sigma**2))*E_cov_k
        return e_g, E_cov_mn


class DensityOfStates(MullikenAnalysis):
    """ A class that calculates different kinds of local and projected
    density of states using the Mulliken charges. """

    def __init__(self, calc):
        MullikenAnalysis.__init__(self, calc)
        self.eigs = self.calc.st.get_eigenvalues()


    def DOS(self, sigma, npts=300):
        """ Return the energy array and corresponding DOS array. """
        eigs = self.eigs
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
        eigs = self.eigs
        e_min, e_max = min(eigs), max(eigs)
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
        eigs = self.eigs
        e_min, e_max = min(eigs), max(eigs)
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
            for k, e_k in enumerate(eigs):
                pdos_k = [nu.exp(-(e-e_k)**2/(2*sigma**2)) for e in e_g]
                for I in indices:
                    q_Ikl = self.mulliken_I_k_l(I,k,li)
                    pdos += nu.array(pdos_k) * q_Ikl
        return e_g, pdos




from ase import Hartree
import numpy as nu
from weakref import proxy
from box.mix import gauss_fct
from box import mix

def get_angular_momenta(l):
    ret = []
    if 's' in l: ret.append(0)
    if 'p' in l: ret.append(1)
    if 'd' in l: ret.append(2)
    return ret


class MullikenAnalysis:
    def __init__(self, calc):
        """
        Class for Mulliken charge analysis.
        
        q_I = sum_k w_k Trace_I Re[ rho(k)*S(k) ]
            = sum_(i in I) Re [ sum_k w_k sum_j rho(k)_ij*S(k)_ji ] 
            = sum_(i in I) Re [ sum_k w_k sum_j rho(k)_ij*S(k)^T_ij ]
            = sum_(i in I) [ sum_k w_k diag(k)_i ],
            = sum_(i in I) diag_i
            
            where diag(k)_i = Re [sum_j rho(k)_ij * S(k)^T_ij] 
            and diag_i = sum_k w_k diag(k)_i
            
        Charge in pieces:
        
        q_(k,a,mu) = wk * f(k,a) * Re [sum_nu wf(k,a,mu)^*wf[k,a,nu]*S(k,mu,nu)]  
                   = wk * f(k,a) * aux(k,a,mu)   
               
            where aux(k,a,mu) = Re [wf(k,a,mu)*sum_nu wf(k,a,nu)*S(k,mu,nu)]     
        
        All the units are, also inside the class, in eV and Angstroms. 
        """
        self.calc = proxy(calc)
        st = self.calc.st
        self.N = self.calc.el.N
        self.nk = self.calc.st.nk
        self.wk = self.calc.st.wk.copy()
        norb = self.calc.st.norb
        self.norb = norb
        
        self.diag = nu.zeros((norb))        
        for k in xrange(self.nk):
            diag_k = ( st.rho[k]*st.S[k].transpose() ).sum(axis=1).real
            self.diag += self.wk[k]*diag_k
            
        self.aux = nu.zeros((self.nk,norb,norb))
        wf = st.wf.copy()
        wfc = st.wf.copy().conjugate()
        for k in xrange(self.nk):
            for a in range(norb):
                self.aux[k,a] = ( wfc[k,a] * nu.dot(wf[k,a],st.S[k].transpose()) ).real 


    def get_rhoa(self,a):
        """ 
        Return the density matrix for eigenstate a, for all k-points.
        """
        rho = nu.zeros_like(self.rhoSk)
        for k,wk in enumerate(self.wk):
            rho[k] = wk*nu.outer(self.calc.st.wf[k,a,:],self.calc.st.wf[k,a,:].conjugate()).real
        return rho


    def trace_I(self,I,matrix):
        """ Return partial trace over atom I's orbitals. """
        ret = 0.0
        I = self.calc.el.orbitals(I,indices=True)
        for i in I:
            ret += matrix[i,i]
        return ret


    def atoms_mulliken(self):
        """ Return Mulliken populations. """
        q = []
        for o1, no in self.calc.el.get_property_lists(['o1','no']):
            q.append( self.diag[o1:o1+no].sum() )
        return nu.array(q)-self.calc.el.get_valences()


    def atom_mulliken(self, I):
        """ 
        Return Mulliken population for atom I.
        
        parameters:
        ===========
        I:        atom index
        """
        orbs = self.calc.el.orbitals(I,indices=True)
        return sum( self.diag[orbs] )


    def basis_mulliken(self, mu):
        """ 
        Return Mulliken population of given basis state.
        
        parameters:
        ===========
        mu:     orbital index (see Elements' methods for indices)
        """
        return self.diag[mu] 
    
    
    def atom_all_angmom_mulliken(self,I):
        """ Return the Mulliken population of atom I from eigenstate a, for all angmom."""
        
        orbs = self.calc.el.orbitals(I,indices=True)
        all = self.diag[orbs]
        pop = nu.zeros((3,))
        pop[0] = all[0]
        if len(all)>1: pop[1] = all[1:4].sum()
        if len(all)>4: pop[2] = all[4:].sum()
        return pop 


    def atom_wf_mulliken(self,I,k,a,wk=True):
        """
        Return Mulliken population for given atom and wavefunction.
               
        parameters:
        ===========
        I:      atom index
        k:      k-vector index
        a:      eigenstate index
        wk:     embed k-point weight in population
        """
        w = 1.0
        if wk: w = self.wk[k]
        orbs = self.calc.el.orbitals(I,indices=True)
        return w*self.aux[k,a,orbs].sum()


    def atom_wf_all_orbital_mulliken(self,I,k,a,wk=True):
        """
        Return orbitals' Mulliken populations for given atom and wavefunction.
        
        parameters:
        ===========
        I:      atom index (returned array size = number of orbitals on I)
        k:      k-vector index 
        a:      eigenstate index
        wk:     embed k-point weight in population
        """
        w = 1.0
        if wk: w = self.wk[k]
        orbs = self.calc.el.orbitals(I,indices=True)
        q = []
        for mu in orbs:
            q.append( w*self.aux[k,a,mu] )
        return nu.array(q)


    def atom_wf_all_angmom_mulliken(self,I,k,a,wk=True):
        """ 
        Return atom's Mulliken populations for all angmom for given wavefunction.
        
        parameters:
        ===========
        I:        atom index
        k:        k-vector index
        a:        eigenstate index
        wk:       embed k-point weight into population
        
        return: array (length 3) containing s,p and d-populations      
        """
        all = self.atom_wf_all_orbital_mulliken(I,k,a,wk)
        pop = nu.zeros((3,))
        pop[0] = all[0]
        if len(all)>1: pop[1] = all[1:4].sum()
        if len(all)>4: pop[2] = all[4:].sum()
        return pop   



class DensityOfStates(MullikenAnalysis):
    def __init__(self, calc):
        """ 
        A class that calculates different kinds of local and projected
        density of states using the Mulliken charges. 
        
        Units also inside this class are in eV
        """
        MullikenAnalysis.__init__(self, calc)

        self.e = self.calc.st.get_eigenvalues()*Hartree
        self.e -= self.calc.get_fermi_level()
        self.eflat = self.e.flatten()
        self.emin, self.emax = self.eflat.min(), self.eflat.max()
        


    def density_of_states(self,broaden=False,width=0.05,window=None,npts=501):
        """
        Return the full density of states.
        
        Sum of states over k-points. Zero is the Fermi-level.
        Spin-degeneracy is NOT counted.
        
        parameters:
        ===========
        broaden:     * If True, return broadened DOS in regular grid
                       in given energy window. 
                     * If False, return energies of all states, followed
                       by their k-point weights. 
        width:       Gaussian broadening (eV)
        window:      energy window around Fermi-energy; 2-tuple (eV)
        npts:        number of data points in output
        """
        mn, mx = self.emin, self.emax
        if window is not None:
            mn, mx = window
            
        x, y = [],[]
        for a in range(self.calc.el.norb):
            x = nu.concatenate( (x,self.e[:,a]) )
            y = nu.concatenate( (y,self.calc.st.wk) )
        x=nu.array(x) 
        y=nu.array(y) 
        if broaden:
            e,dos = mix.broaden(x, y, width=width, N=npts, a=mn, b=mx)
        else:
            e,dos = x,y
        return e,dos


    def local_density_of_states(self,projected=False,width=0.05,window=None,npts=501):
        """
        Return state density for all atoms as a function of energy.
        
        parameters:
        ===========
        projected: return local density of states projected for 
                   angular momenta 0,1 and 2 (s,p and d)
                   ( sum of pldos over angular momenta = ldos ) 
        width:     energy broadening (in eV)
        window:    energy window around Fermi-energy; 2-tuple (eV)
        npts:      number of grid points for energy
        
        return:    projected==False:
                        energy grid, ldos[atom,grid]
                   projected==True:
                        energy grid, 
                        ldos[atom, grid],
                        pldos[atom, angmom, grid]
        """
        mn, mx = self.emin, self.emax
        if window is not None:
            mn, mx = window
                   
        # only states within window 
        kl,al,el = [],[],[]        
        for k in range(self.nk):
            for a in range(self.norb):
                if mn<=self.e[k,a]<=mx:
                    kl.append(k)
                    al.append(a)
                    el.append(self.e[k,a]) 
                
        ldos = nu.zeros((self.N,npts))
        if projected:
            pldos = nu.zeros((self.N,3,npts))
        for i in range(self.N):
            q = [ self.atom_wf_mulliken(i,k,a,wk=True) for k,a in zip(kl,al) ]
            egrid, ldos[i,:] = mix.broaden( el,q,width=width,N=npts,a=mn,b=mx )
            if projected:
                q = nu.array( [self.atom_wf_all_angmom_mulliken(i,k,a,wk=True) for k,a in zip(kl,al)] )
                for l in range(3):
                    egrid, pldos[i,l] = mix.broaden( el,q[:,l],width=width,N=npts,a=mn,b=mx )
            
        if projected:
            assert nu.all( abs(ldos-pldos.sum(axis=1))<1E-6 )
            return egrid, ldos, pldos
        else:       
            return egrid,ldos


    def projected_local_density_of_states(self,width=0.05,window=None,npts=501):
        """
        Return state density for all atoms as a function of energy and angmom.
        
        parameters:
        ===========
        width:     energy broadening (in eV)
        window:    energy window around Fermi-energy; 2-tuple (eV)
        npts:      number of grid points for energy
        
        return:    energy grid, ldos[atom-index, angmom, energy grid-index]
        """
        mn, mx = self.emin, self.emax
        if window is not None:
            mn, mx = window
                   
        # only states within window 
        kl,al,el = [],[],[]        
        for k in range(self.nk):
            for a in range(self.norb):
                if mn<=self.e[k,a]<=mx:
                    kl.append(k)
                    al.append(a)
                    el.append(self.e[k,a]) 
                
        pldos = nu.zeros((self.N,3,npts))
        for i in range(self.N):
            q = nu.array( [self.atom_wf_all_angmom_mulliken(i,k,a,wk=True) for k,a in zip(kl,al)] )
            for l in range(3):
                egrid, pldos[i,l] = mix.broaden( el,q[:,l],width=width,N=npts,a=mn,b=mx )
        return egrid,pldos

        


class MullikenBondAnalysis(MullikenAnalysis):    
    
    def __init__(self, calc):
        """ 
        Class for bonding analysis using Mulliken charges. 
        """
        raise NotImplementedError('needs re-writing for k-points')
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
            rho_k_mn = self.calc.st.wf[0,k,m]*self.calc.st.wf[0,k,n].conjugate()
            if occupations:
                rho_k_mn *= f[k]
            E_cov_k = rho_k_mn*mat_nm
            add = self.gauss_fct(e_g, e_k, sigma)*E_cov_k
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
            add = self.gauss_fct(e_g, e_k, sigma)
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
        dq_I = self.atom_mulliken(I) - self.calc.el.get_valences()[I]
        return 0.5*gamma_II*dq_I**2 + self.E_prom_I(I)


    def E_prom_I(self, I):
        """ Return the promotion energy of atom I. """
        orb_indices = self.calc.el.orbitals(I, indices=True)
        ret = 0.0
        for m in orb_indices:
            q_mu = self.basis_mulliken(m)
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
        dq_I = self.atom_mulliken(I) - self.calc.el.get_valences()[I]
        dq_J = self.atom_mulliken(J) - self.calc.el.get_valences()[J]

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



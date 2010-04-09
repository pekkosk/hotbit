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
        self.st = st
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


    def trace_I(self,I,matrix):
        """ Return partial trace over atom I's orbitals. """
        ret = 0.0
        I = self.calc.el.orbitals(I,indices=True)
        for i in I:
            ret += matrix[i,i]
        return ret


    def get_atoms_mulliken(self):
        """ Return Mulliken populations. """
        q = []
        for o1, no in self.calc.el.get_property_lists(['o1','no']):
            q.append( self.diag[o1:o1+no].sum() )
        return nu.array(q)-self.calc.el.get_valences()


    def get_atom_mulliken(self, I):
        """ 
        Return Mulliken population for atom I.
        
        parameters:
        ===========
        I:        atom index
        """
        orbs = self.calc.el.orbitals(I,indices=True)
        return sum( self.diag[orbs] )


    def get_basis_mulliken(self, mu):
        """ 
        Return Mulliken population of given basis state.
        
        parameters:
        ===========
        mu:     orbital index (see Elements' methods for indices)
        """
        return self.diag[mu] 
    
    
    def get_atom_all_angmom_mulliken(self,I):
        """ Return the Mulliken population of atom I from eigenstate a, for all angmom."""
        
        orbs = self.calc.el.orbitals(I,indices=True)
        all = self.diag[orbs]
        pop = nu.zeros((3,))
        pop[0] = all[0]
        if len(all)>1: pop[1] = all[1:4].sum()
        if len(all)>4: pop[2] = all[4:].sum()
        return pop 


    def get_atom_wf_mulliken(self,I,k,a,wk=True):
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


    def get_atom_wf_all_orbital_mulliken(self,I,k,a,wk=True):
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


    def get_atom_wf_all_angmom_mulliken(self,I,k,a,wk=True):
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
        all = self.get_atom_wf_all_orbital_mulliken(I,k,a,wk)
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
        


    def get_density_of_states(self,broaden=False,width=0.05,window=None,npts=501):
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


    def get_local_density_of_states(self,projected=False,width=0.05,window=None,npts=501):
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
            q = [ self.get_atom_wf_mulliken(i,k,a,wk=True) for k,a in zip(kl,al) ]
            egrid, ldos[i,:] = mix.broaden( el,q,width=width,N=npts,a=mn,b=mx )
            if projected:
                q = nu.array( [self.get_atom_wf_all_angmom_mulliken(i,k,a,wk=True) for k,a in zip(kl,al)] )
                for l in range(3):
                    egrid, pldos[i,l] = mix.broaden( el,q[:,l],width=width,N=npts,a=mn,b=mx )
            
        if projected:
            assert nu.all( abs(ldos-pldos.sum(axis=1))<1E-6 )
            return egrid, ldos, pldos
        else:       
            return egrid,ldos

       


class MullikenBondAnalysis(MullikenAnalysis):    
    
    def __init__(self, calc):
        """ 
        Class for bonding analysis using Mulliken charges. 
        """
        #raise NotImplementedError('needs re-writing for k-points')
        MullikenAnalysis.__init__(self, calc)
        rhoSk = []
        for k in range(self.nk):
            rhoSk.append( nu.dot(self.st.rho[k],self.st.S[k]) ) 
        self.rhoSk = nu.array(rhoSk)
        self.SCC = self.calc.get('SCC')
        
        
    def get_mayer_bond_order(self,i,j):
        """
        Return Mayer bond-order between two atoms.
        
        Warning: appears to work only with periodic systems
        where orbitals have no overlap with own images.
        
        parameters:
        ===========
        I:        first atom index
        J:        second atom index
        """
        assert type(i)==int and type(j)==int
        orbi = self.calc.el.orbitals(i, indices=True)
        orbj = self.calc.el.orbitals(j, indices=True)
        
        M = 0.0
        for k in xrange(self.nk):
            for mu in orbi:
                for nu in orbj:
                    M += self.wk[k] * self.rhoSk[k,mu,nu]*self.rhoSk[k,nu,mu]
        assert abs(M.imag)<1E-12
        return M.real


    def get_atom_energy(self, I):
        """ 
        Return the absolute atom energy (in eV).
        
        parameters:
        ===========
        I:         atom index
        """
        if self.SCC:
            coul = 0.5*self.calc.st.es.G[I,I]*self.st.dq[I]**2
        else:
            coul = 0.0         
        rep = self.calc.rep.get_pair_repulsive_energy(I,I) #self-repulsion for pbc
        return coul*Hartree + self.get_promotion_energy(I) + rep*Hartree


    def get_promotion_energy(self, I):
        """ 
        Return atom's promotion energy (in eV). 
        
        energy = sum_k w_k sum_(m,n in I) rho_mn(k) H^0_nm(k) - E_free(I)
        
        parameters:
        ===========
        I:         atom index
        """
        orbi = self.calc.el.orbitals(I, indices=True)
        e = 0.0
        for k in range(self.nk):
            for m1 in orbi:
                for m2 in orbi:
                    e += self.wk[k]*self.st.rho[k,m1,m2]*self.st.H0[k,m2,m1]
        assert abs(e.imag)<1E-12
        e = e - self.calc.el.elements[self.calc.el.symbols[I]].get_free_atom_energy()
        return e.real * Hartree


    def get_bond_energy(self,i,j):
        """ 
        Return the absolute bond energy between atoms (in eV). 
        
        parameters:
        ===========
        i,j:     atom indices
        """
        rep = self.calc.rep.get_pair_repulsive_energy(i,j)
        if self.SCC:
            coul = self.st.es.G[i,j]*self.st.dq[i]*self.st.dq[j]
        else:
            coul = 0.0
                
        orbi = self.calc.el.orbitals(i,indices=True)
        orbj = self.calc.el.orbitals(j,indices=True)
        ebs = 0.0
        for k in range(self.nk):
            for mu in orbi:
                for nu in orbj:
                    ebs += self.wk[k]* 2*(self.st.rho[k,mu,nu]*self.st.H0[k,nu,mu] ).real
        return (rep + coul + ebs)*Hartree


    def get_atom_and_bond_energy(self, i):
        """
        Return given atom's contribution to cohesion (in eV).
        
        parameters:
        ===========
        i:    atom index
        """
        ea = self.get_atom_energy(i)
        eb = 0.0
        for j in range(self.N):
            if j==i: continue
            eb += 0.5 * self.get_bond_energy(i,j)
        return ea + eb


    def get_covalent_energy(self):
        """ Returns the covalent bond energy of the whole system. """
        E_bs = self.calc.st.band_structure_energy()*Hartree
        E = nu.sum(self.rho*self.bar_epsilon*self.S)
        return E_bs - E


    def get_E_cov_m_n(self, m, n, sigma, npts=500, e_min=None, e_max=None, occupations=False):
        # name-> orbital_pair_covalent_energy
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
        # name: atom_pair_covalent_energy
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
        # angmom_pair_covalent_energy
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
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
        


    def get_density_of_states(self,broaden=False,projected=False,occu=False,width=0.05,window=None,npts=501):
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
        projected:   project DOS for angular momenta 
        occu:        for not broadened case, return also state occupations
        width:       Gaussian broadening (eV)
        window:      energy window around Fermi-energy; 2-tuple (eV)
        npts:        number of data points in output
        
        return:      * if projected: e[:],dos[:],pdos[l,:] (angmom l=0,1,2)
                     * if not projected: e[:],dos[:]
                       * if broaden: e[:] is on regular grid, otherwise e[:] are
                         eigenvalues and dos[...] corresponding weights
                     * if occu: e[:],dos[:],occu[:] 
                          
        """
        if broaden and occu or projected and occu:
            raise AssertionError('Occupation numbers cannot be returned with broadened DOS.')
        if projected:
            e,dos,pdos = self.get_local_density_of_states(True,width,window,npts)
            return e,dos.sum(axis=0),pdos.sum(axis=0)
            
        mn, mx = self.emin, self.emax
        if window is not None:
            mn, mx = window
            
        x, y, f = [],[],[]
        for a in range(self.calc.el.norb):
            x = nu.concatenate( (x,self.e[:,a]) )
            y = nu.concatenate( (y,self.calc.st.wk) )
            f = nu.concatenate( (f,self.calc.st.f[:,a]) )
        x=nu.array(x) 
        y=nu.array(y) 
        f=nu.array(f)
        if broaden:
            e,dos = mix.broaden(x, y, width=width, N=npts, a=mn, b=mx)
        else:
            e,dos = x,y
        if occu:
            return e,dos,f
        else:
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
            return egrid, ldos

       


class MullikenBondAnalysis(MullikenAnalysis):    
    
    def __init__(self, calc):
        """ 
        Class for bonding analysis using Mulliken charges. 
        """
        #raise NotImplementedError('needs re-writing for k-points')
        MullikenAnalysis.__init__(self, calc)
        rhoSk = []
        HS = []
        rhoM = []
        n=self.st.norb
        aux = nu.diag(self.st.H0[0]).repeat(n).reshape((n,n))
        #epsilon-bar matrix of Bornsen et al J.Phys.:Cond.mat 11, L287 (1999)
        eps = 0.5*(aux+aux.transpose())   
        for k in range(self.nk):
            rhoSk.append( nu.dot(self.st.rho[k],self.st.S[k]) )
            HS.append( self.st.H0[k] - self.st.S[k]*eps )
            rhoM.append( self.st.rho[k]*((self.st.H0[k]-self.st.S[k]*eps).transpose()) ) 
        self.rhoSk = nu.array(rhoSk)
        self.rhoM = nu.array(rhoM)
        self.HS = nu.array(HS)
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
        
        Warning: bonding & atom energy analysis less clear for
        systems where orbitals overlap with own periodic images.
        
        parameters:
        ===========
        I:         atom index
        """
        if self.SCC:
            coul = 0.5*self.calc.st.es.G[I,I]*self.st.dq[I]**2
        else:
            coul = 0.0         
        
        elm = self.calc.el.elements[self.calc.el.symbols[I]]
        o1, no = self.calc.el.get_property_lists(['o1','no'])[I]
        
        ecorr = 0.0
        for mu in range(o1,o1+no):
            emu = self.calc.el.orbitals(mu,basis=True)['energy']
            ecorr += self.calc.el.get_free_population(mu) * (self.st.H0[0,mu,mu]-emu)
            #ecorr += self.get_basis_mulliken(mu)*self.st.H0[0,mu,mu] - self.calc.el.get_free_population(mu)*emu
        
        eorb = 0.0
        for k,wk in enumerate(self.wk):
            eorb += wk*nu.sum( self.rhoM[k,o1:o1+no,o1:o1+no] )
        #eorb = sum( [self.wk[k]*nu.sum(self.rhoM[k,o1:o1+no,o1:o1+no]) for k in range(self.nk)] )
        
        erep = self.calc.rep.get_pair_repulsive_energy(I,I) #self-repulsion for pbc
        A = (coul + erep + ecorr + eorb )*Hartree + self.get_promotion_energy(I) 
        assert abs(A.imag)<1E-12
        return A.real


    def get_promotion_energy(self, I):
        """ 
        Return atom's promotion energy (in eV). 
        
        Defined as:
            E_prom,I = sum_(mu in I) [q_(mu) - q_(mu)^0] H_(mu,mu)(k=0)
        
        Warning: bonding & atom energy analysis less clear for
        systems where orbitals overlap with own periodic images.
        
        parameters:
        ===========
        I:         atom index
        """
        
        #FIXME H0[0,:,:] is often not the gamma-point!!!
        e = 0.0
        for mu in self.calc.el.orbitals(I, indices=True):
            q = self.get_basis_mulliken(mu)
            q0 = self.calc.el.get_free_population(mu)
            e += (q-q0)*self.st.H0[0,mu,mu]
            
        assert abs(e.imag)<1E-12
        return e.real * Hartree


    def get_bond_energy(self,i,j):
        """ 
        Return the absolute bond energy between atoms (in eV).
        
        Warning: bonding & atom energy analysis less clear for
        systems where orbitals overlap with own periodic images. 
        
        parameters:
        ===========
        i,j:     atom indices
        """
        rep = self.calc.rep.get_pair_repulsive_energy(i,j)
        if self.SCC:
            coul = self.st.es.G[i,j]*self.st.dq[i]*self.st.dq[j]
        else:
            coul = 0.0
                
        o1i, noi = self.calc.el.get_property_lists(['o1','no'])[i]
        o1j, noj = self.calc.el.get_property_lists(['o1','no'])[j]
        ebs = 0.0
        for k,wk in enumerate(self.wk):
            ebs += 2*wk*nu.sum( self.rhoM[k,o1i:o1i+noi,o1j:o1j+noj].real )
        return (rep + coul + ebs) * Hartree


    def get_atom_and_bond_energy(self, i):
        """
        Return given atom's contribution to cohesion (in eV).
        
        Warning: bonding & atom energy analysis less clear for
        systems where orbitals overlap with own periodic images.
        
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

    
    def get_covalent_energy(self,mode='default',i=None,j=None,width=None,window=None,npts=501):
        """
        Return covalent bond energies in different modes. (eV)
        
        ecov is described in 
        Bornsen, Meyer, Grotheer, Fahnle, J. Phys.:Condens. Matter 11, L287 (1999) and
        Koskinen, Makinen Comput. Mat. Sci. 47, 237 (2009)
        
        
        
        parameters:
        ===========
        mode:    'default' total covalent energy
                 'orbitals' covalent energy for orbital pairs
                 'atoms' covalent energy for atom pairs
                 'angmom' covalent energy for angular momentum components
        i,j:     atom or orbital indices, or angular momentum pairs
        width:   * energy broadening (in eV) for ecov
                 * if None, return energy eigenvalues and corresponding 
                   covalent energies in arrays, directly
        window:  energy window (in eV wrt Fermi-level) for broadened ecov
        npts:    number of points in energy grid (only with broadening) 
    
        return:
        =======
        x,y:     * if width==None, x is list of energy eigenvalues (including k-points)
                   and y covalent energies of those eigenstates
                 * if width!=None, x is energy grid for ecov.
                 * energies (both energy grid and ecov) are in eV.
         
        Note: energies are always shifted so that Fermi-level is at zero. 
              Occupations are not otherwise take into account (while k-point weights are)
        """
        eps = 1E-6
        x, y = [], []
        wf = self.st.wf
        energy = self.st.e - self.st.occu.get_mu()
        if window==None:
            mn,mx = energy.flatten().min()-eps, energy.flatten().max()+eps
        else:
            mn,mx = window[0]/Hartree, window[1]/Hartree
        
        if mode=='angmom':
            lorbs = [[],[],[]]
            for m,orb in enumerate(self.calc.el.orbitals()):
                lorbs[orb['angmom']].append(m)
            oi, oj = nu.array(lorbs[i]), nu.array(lorbs[j])
        elif mode=='atoms':
            o1i, noi = self.calc.el.get_property_lists(['o1','no'])[i]
            o1j, noj = self.calc.el.get_property_lists(['o1','no'])[j]
            
        for k,wk in enumerate(self.wk):
            for a in range(self.norb):
                if not mn<=energy[k,a]<=mx:
                    continue 
                x.append( energy[k,a] )
                if mode == 'default':
                    e = 0.0
                    for m in range(self.norb):
                        e += wk*nu.sum( (wf[k,a,m]*wf[k,a,:].conj()*self.HS[k,:,m]) )
                    y.append(e)
                elif mode == 'orbitals':
                    if i!=j:
                        y.append( wk*2*(wf[k,a,i]*wf[k,a,j].conj()*self.HS[k,j,i]).real )
                    else:
                        y.append( wk*wf[k,a,i]*wf[k,a,j].conj()*self.HS[k,j,i])
                elif mode == 'atoms':
                    e = 0.0
                    if i!=j:
                        for m in range(o1i,o1i+noi):
                            for n in range(o1j,o1j+noj):
                                e += wk*2*(wf[k,a,m]*wf[k,a,n].conj()*self.HS[k,n,m]).real
                    else:
                        for m in range(o1i,o1i+noi):
                            for n in range(o1j,o1j+noj):
                                e += wk*(wf[k,a,m]*wf[k,a,n].conj()*self.HS[k,n,m])
                    y.append(e)
                elif mode == 'angmom':
                    e = 0.0                                
                    for m in lorbs[i]:                           
                        e += wk*nu.sum( wf[k,a,m]*wf[k,a,oj].conj()*self.HS[k,oj,m] )
                    if i!=j:
                        e += e.conj()
                    y.append(e)
                else:
                    raise NotImplementedError('Unknown covalent energy mode "%s".' %mode) 
                    
        x,y = nu.array(x), nu.array(y)
        assert nu.all( abs(y.imag)<1E-12 )
        y=y.real
        if width==None:
            return x * Hartree, y * Hartree
        else:
            e,ecov = mix.broaden(x, y, width=width/Hartree, N=npts, a=mn, b=mx)
            return e * Hartree, ecov * Hartree
                   
                

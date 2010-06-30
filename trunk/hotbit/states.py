# Please see the accompanying LICENSE file for further information.

from solver import Solver
from weakref import proxy
from electrostatics import Electrostatics
from occupations import Occupations
import numpy as np
from box import mix
from auxil import k_to_kappa_points
from box.mix import divisors
pi=np.pi


#
# Constructing the density matrix, and the energy
# weighted density matrix from eigenvectors, eigenvalues
# and occupations.
#

def compute_rho(wf, occ):
    nk, n, m = wf.shape
    nk2, k = occ.shape

    assert n == m
    assert n == k
    assert nk == nk2

    rho = np.zeros(wf.shape, dtype=wf.dtype)
    for k in range(nk):
        rho[k,:,:] = np.dot(
            wf[k,:,:].transpose(),
            occ[k,:].reshape(n,-1) * wf[k,:,:].conj()
            )
    return rho


def compute_rhoe(wf, occ, e):
    nk, n, m = wf.shape
    nk2, k = occ.shape

    assert n == m
    assert n == k
    assert nk == nk2

    rhoe = np.zeros(wf.shape, dtype=wf.dtype)
    for k in range(nk):
        rhoe[k,:,:] = np.dot(
            wf[k,:,:].transpose(),
            e[k,:].reshape(n,-1)*occ[k,:].reshape(n,-1) * wf[k,:,:].conj()
            )
    return rhoe


#
# States class
#

class States:

    def __init__(self,calc):
        self.es = Electrostatics(calc,
            charge_density = calc.get('charge_density'),
            solver = calc.get('coulomb_solver'))
        self.solver=Solver(calc)
        self.calc=proxy(calc)
        self.nat=len(calc.el)
        self.norb=calc.el.get_nr_orbitals()
        self.prev_dq=[None,None]
        self.count=0
        self.SCC=calc.get('SCC')
        self.rho=None
        self.rhoe0=None
        self.nk=None
        
       
    def setup_k_sampling(self,kpts,physical=True,rs='kappa'):
        '''
        Setup the k-point sampling and their weights.
        
        @param kpts: 3-tuple: number of k-points in different directions
                     list of 3-tuples: k-points given explicitly
        @param physical: No meaning for infinite periodicities. For physically
                     periodic systems only certain number of k-points are allowed.
                     (like wedge of angle 2*pi/N, only number of k-points that 
                     divides N, is physically allowed). If physical=False,
                     allow interpolation of this k-sampling.
        rs:          use 'kappa'- or 'k'-point sampling
        '''
        if kpts!=(1,1,1) and self.calc.get('width')<1E-10:
            raise AssertionError('With k-point sampling width must be>0!')
            
        M = self.calc.el.get_number_of_transformations()
        if isinstance(kpts,tuple):
            table = self.calc.el.atoms.container.table
            # set up equal-weighted and spaced k-point mesh
            if 0 in kpts:
                raise AssertionError('Each direction must have at least one k-point! (Gamma-point)')
            
            kl=[]
            for i in range(3):
                if M[i]==np.Inf:
                    # arbitrary sampling is allowed
                    spacing = 2*pi/kpts[i]
                    kl.append( np.linspace(-pi+spacing/2,pi-spacing/2,kpts[i]) )
                else:
                    
                    # TODO: choose the closest number of k-points if
                    # sampling should be physical
                    # first calculate possible numer of k-points, then
                    # select the ones closest to the desired one
                    # nks = M/divisors(M) 
                    # nk1 = nks[abs(nks-nk).argmin()]                  
                    
                    # discrete, well-defined sampling; any k-point is not allowed 
                    if kpts[i] not in mix.divisors(M[i]) and physical:
                        print 'Allowed k-points for direction',i,'are',mix.divisors(M[i])
                        raise Warning('Non-physical k-point sampling! ')
                    else:
                        kl.append( np.linspace(0,2*pi-2*pi/kpts[i],kpts[i]) )
                        
                        
                
            k=[]    
            wk=[]
            nk0 = np.prod(kpts)
            for a in range(kpts[0]):
                for b in range(kpts[1]):
                    for c in range(kpts[2]):
                        newk = np.array([kl[0][a], kl[1][b], kl[2][c]])
                        one_equivalent = False
                        for i in range(3):
                            if 'equivalent' in table[i]:
                                if one_equivalent:
                                    raise NotImplementedError('Surprising type of new symmetry; reconsider implementation.')
                                one_equivalent = True
                                n = table[i]['equivalent'] # symmetry op i is equivalent to tuple n
                                assert n[i]==0
                                #for j in range(3):
                                #    if n[j]!=0: break
                                #assert i!=j
                                #assert abs(n[j])>0
                                #newk[i] = newk[i] + 1.0*n[j]/M[i]*newk[j]
                                newk[i] = newk[i] + 1.0*np.dot(n,newk)/M[i] #n[j]/M[i]*newk[j]
                                
                                
                        inv_exists = False
                        # if newk's inverse exists, increase its weight by default
                        for ik, oldk in enumerate(k):
                            if np.linalg.norm(oldk+newk)<1E-10: 
                                inv_exists = True
                                wk[ik]+=1.0/nk0
                        # newk's inverse does not exist; make new k-point
                        if not inv_exists:
                            k.append( newk )
                            wk.append( 1.0/nk0 ) 
            nk=len(k)            
            k=np.array(k)
            wk=np.array(wk)
        else:
            # work with a given set of k-points
            nk=len(kpts)
            if rs=='k':
                k = k_to_kappa_points(kpts,self.calc.el.atoms)
            else:
                k=np.array(kpts)
            wk=np.ones(nk)/nk
            kl=None
            
            
        # now sampling is set up. Check the consistency.
        pbc = self.calc.el.get_pbc()
        for i in range(3):
            for kp in k:
                if kp[i]>1E-10 and not pbc[i]:
                    raise AssertionError('Do not set (non-zero) k-points in non-periodic direction!')
        
        return nk, k, kl, wk
        

    def guess_dq(self):
        n=len(self.calc.el)
        if not self.SCC:
            return np.zeros((n,))
        if self.count==0:
            return np.zeros((n,)) - float(self.calc.get('charge'))/n
        elif self.count==1:    # use previous charges
            return self.prev_dq[0]
        else:                  # linear extrapolation
            return self.prev_dq[0] + (self.prev_dq[0]-self.prev_dq[1])


    def solve(self):
        if self.nk==None:
            physical = self.calc.get('physical_k')
            self.nk, self.k, self.kl, self.wk = self.setup_k_sampling( self.calc.get('kpts'),physical=physical,rs=self.calc.get('rs') )
            width=self.calc.get('width')
            self.occu = Occupations(self.calc.el.get_number_of_electrons(),width,self.wk)
        self.calc.start_timing('solve')
        
        # TODO: enable fixed dq-calculations in SCC (for band-structures)
        dq=self.guess_dq()
        if self.calc.dry_run:
            self.calc.memory_estimate()
        self.H0, self.S, self.dH0, self.dS = self.calc.ia.get_matrices()
#        try:
        if True:
            if self.SCC:
                self.es.construct_Gamma_matrix(self.calc.el.atoms)
            self.e, self.wf = self.solver.get_states(self.calc,dq,self.H0,self.S)
            
            self.check_mulliken_charges()
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
            print>> self.calc.get_output(), "Mulliken populations: min=%0.3f, max=%0.3f" % (np.min(dQ),np.max(dQ))
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
        self.rho = compute_rho(self.wf,self.f)
        self.calc.stop_timing('rho')
        if self.SCC:
            self.dq = self.mulliken()
            self.es.set_dq(self.dq)
        self.calc.stop_timing('update')


    def large_update(self):
        """ Update stuff from eigenstates needed later for forces etc. """
        self.calc.start_timing('final update')
        self.Swf = None
        self.dH = np.zeros_like(self.dH0)
        if self.SCC:
            self.prev_dq=[self.dq, self.prev_dq[0]]
            # TODO: do k-sum with numpy
            for ik in range(self.nk):
                for a in range(3):
                    self.dH[ik,:,:,a] = self.dH0[ik,:,:,a] + self.es.get_h1()*self.dS[ik,:,:,a]
        else:
            self.dH = self.dH0

        # density matrix weighted by eigenenergies
        self.rhoe = compute_rhoe(self.wf,self.f,self.e)
        self.calc.stop_timing('final update')


    def get_dq(self):
        return self.dq.copy()


    def get_eigenvalues(self):
        return self.e.copy()


    def get_band_energies(self, kpts):
        """
        Return the eigenvalue spectrum for a set of explicitly given k-points.
        """
        H0, S, dH0, dS = self.calc.ia.get_matrices(kpts)
        e, wf = self.solver.get_eigenvalues_and_wavefunctions(H0, S)

        return e


    def get_occupations(self):
        return self.f.copy()


    def get_homo(self,occu=0.99):
        """ Return highest *largely* occupied orbital (>occu)
        
        0<=occu<=2. Recommended use only for molecules.
        """
        for i in range(self.norb)[::-1]:
            if np.any( self.f[:,i]>occu ): return i


    def get_lumo(self,occu=1.01):
        """ Return lowest *largely* unuccopied orbital (<occu)
        
        0<=occu<=2. Recommended use only for molecules.
        """
        for i in range(self.norb):
            if np.any( self.f[:,i]<occu ): return i


    def mulliken(self):
        '''
        Return excess Mulliken populations dq = dq(total)-dq0
        
        dq_I = sum_k w_k Trace_I Re[ rho(k)*S(k) ]
             = sum_(i in I) Re [ sum_k w_k sum_j rho(k)_ij*S(k)_ji ] 
             = sum_(i in I) Re [ sum_k w_k sum_j rho(k)_ij*S(k)^T_ij ]
             = sum_(i in I) [ sum_k w_k diag(k)_i ],
             = sum_(i in I) diag_i
            
               where diag(k)_i = Re [sum_j rho(k)_ij * S(k)^T_ij] 
               and diag_i = sum_k w_k diag(k)_i 
        '''
        diag = np.zeros((self.norb))
        for ik in xrange(self.nk):
            diag_k = ( self.rho[ik]*self.S[ik].transpose() ).sum(axis=1).real
            diag = diag + self.wk[ik] * diag_k
        q=[]
        for o1, no in self.calc.el.get_property_lists(['o1','no']):
            q.append( diag[o1:o1+no].sum() )
        return np.array(q)-self.calc.el.get_valences()


    def get_band_structure_energy(self):
        '''
        Return band structure energy.
        
        ebs = sum_k w_k ( sum_ij rho_ij * H0_ji )
            = sum_k w_k ( sum_i [sum_j rho_ij * H0^T_ij] )
        '''
        self.calc.start_timing('e_bs')
        ebs = 0.0
        for ik in xrange(self.nk):
            diagonal = ( self.rho[ik]*self.H0[ik].transpose() ).sum(axis=1)
            ebs += self.wk[ik] * diagonal.sum()
        assert ebs.imag<1E-13
        self.calc.stop_timing('e_bs')
        return ebs.real 


    def get_band_structure_forces(self):
        '''
        Return band structure forces.
        
        F_I = - sum_k w_k Trace_I [ dH(k)*rho(k) - dS(k)*rhoe(k) + c.c ]
            = - sum_k w_k sum_(i in I) diag_i(k) + c.c.,
            
                where diag_i(k) = [dH(k)*rho(k) - dS(k)*rhoe(k)]_ii
                                = sum_j [dH(k)_ij*rho(k)_ji - dS(k)_ij*rhoe(k)_ji]
                                = sum_j [dH(k)_ij*rho(k)^T_ij - dS(k)_ij*rhoe(k)^T_ij]
        '''
        self.calc.start_timing('f_bs')       
        diag = np.zeros((self.norb,3),complex)
        
        for a in range(3):
            for ik in range(self.nk):
                diag_k = ( self.dH[ik,:,:,a]*self.rho[ik].transpose()  \
                         - self.dS[ik,:,:,a]*self.rhoe[ik].transpose() ).sum(axis=1)
                diag[:,a] = diag[:,a] - self.wk[ik] * diag_k
            
        f=[]            
        for o1, no in self.calc.el.get_property_lists(['o1','no']):
            f.append( 2*diag[o1:o1+no,:].sum(axis=0).real )
            
        self.calc.stop_timing('f_bs')
        return f


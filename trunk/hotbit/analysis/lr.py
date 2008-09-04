import numpy as nu
from numpy.linalg.linalg import eigh 
from ase.units import Hartree
dot=nu.dot
sqrt=nu.sqrt
import sys
from box.timing import Timer

class LinearResponse:
    """ 
    Calculate linear response in the spirit of LR-TD-DFT.
    
    For details, see Niehaus et.al. Phys. Rev. B 63, 085108 (2001)
    """
    def __init__(self,calc,energy_cut=1.0,timing=False,out=None):
        """ Construct the object.
        
        parameters:
        -----------
        calc: calculator object
        energy_cut: max energy for particle-hole excitations
        timing: output timing summary after calculation
        out: output object (file name or object)
        """
        self.calc=calc
        self.st=calc.st
        self.el=calc.el
        self.es=calc.es
        self.energy_cut=energy_cut
        self.noc=self.st.get_hoc()+1 #number of occupied states (not index)
        self.nel=self.el.get_number_of_electrons()
        self.norb=self.el.get_nr_orbitals()
        self.e=self.st.get_eigenvalues()
        self.N=len(self.el)
        
        if abs(nu.mod(self.nel,2))>1E-2:
            raise RuntimeError('Linear response only for closed shell systems! (even number of electrons)')
        if abs(self.nel-2*self.noc)>1E-2:
            print 'Number of electrons:',self.nel
            print '2*Number of occupied states:',2*self.noc
            raise RuntimeError('Number of electrons!=2*number of occupied orbitals. Decrease electronic temperature?')
        if out is None:
            self.out=sys.stdout
        else:
            self.out=open(out,'a')
        self.timer=Timer('Linear Response',self.out)   
        self.timing=timing         
            
        
    def get_linear_response(self):
        """ Get linear response spectrum in eV. """
        return self.omega*Hartree,self.F        
        
        
    def run(self):
        print>>self.out, '\nLR for %s (charge %.2f). ' %(self.el.atoms.get_name(),self.calc.get_charge()),
        # select electron-hole excitations (i occupied, j not occupied)
        de=[]
        particle_holes=[]
        self.timer.start('setup ph pairs')
        for i in range(self.noc):
            for j in range(self.noc,self.norb):
                energy=self.e[j]-self.e[i]
                if energy<self.energy_cut:
                    particle_holes.append([i,j])
                    de.append(energy)
        self.timer.stop('setup ph pairs')
        
        # setup the matrix (gamma-approximation) and diagonalize
        self.timer.start('setup matrix')
        dim=len(de)
        print>>self.out, 'Dimension %i. ' %dim,
        if not 0<dim<50000:
            raise RuntimeError('Coupling matrix too large or small (%i)' %dim)
        r=self.el.get_positions()            
        transfer_q=nu.array([self.st.mulliken_transfer(ph[0],ph[1]) for ph in particle_holes])
        rv=nu.array([dot(tq,r) for tq in transfer_q])
        gamma=nu.zeros((self.N,self.N))
        for i in range(self.N):
            for j in range(self.N):
                gamma[i,j]=self.es.gamma(i,j)
        
        matrix=nu.zeros((dim,dim))
        gamma_tq=nu.zeros((dim,self.N))
        for k in range(dim):
            gamma_tq[k,:]=dot(gamma,transfer_q[k,:])
            
        for k1,ph1 in enumerate(particle_holes):
            matrix[k1,k1]=de[k1]**2
            for k2,ph2 in enumerate(particle_holes):
                coupling=dot(transfer_q[k1,:],gamma_tq[k2,:])
                matrix[k1,k2]+=2*sqrt(de[k1]*de[k2])*coupling
        self.timer.stop('setup matrix')                            
                                                                        
        print>>self.out, 'coupling matrix constructed. ',
        self.out.flush()         
        self.timer.start('diagonalize')                                                         
        omega2,eigv=eigh(matrix)
        self.timer.stop('diagonalize')
        print>>self.out, 'Matrix diagonalized.',
        self.out.flush()
        assert all(omega2>1E-9)           
        omega=sqrt(omega2)
        
        # calculate oscillator strengths
        F=[]
        self.timer.start('oscillator strengths')
        for ex in range(dim):
            v=[]    
            for i in range(3):
                v.append( sum( rv[:,i]*sqrt(de[:])*eigv[:,ex])/sqrt(omega[ex])*2 )
            F.append( omega[ex]*dot(v,v)*2.0/3 )
        self.omega=omega
        self.F=F   
        self.eigv=eigv
        self.dim=dim
        self.particle_holes=particle_holes
        self.timer.stop('oscillator strengths')
        if self.timing:
            self.timer.summary()
            
    def print_info(self):
        print 'eV..'
        for ex in range(self.dim):
            print self.omega[ex]*Hartree, 
            order=nu.argsort(abs(self.eigv[:,ex]))[::-1]
            for ph in order[:4]:
                i,j=self.particle_holes[ph]
                print '%3i-%-3i:%-10.3f' %(i,j,self.eigv[ph,ex]**2),
            print                 
                
                
                
                
                    
        
        
        
        
        


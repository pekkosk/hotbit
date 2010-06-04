import numpy as nu
from numpy.linalg.linalg import eigh 
from box import mix
from ase.units import Hartree
dot=nu.dot
sqrt=nu.sqrt
import sys
from box.timing import Timer

class LinearResponse:
    """ 
    
    """
    def __init__(self,calc,energy_cut=10.0,timing=False,txt=None):
        """
        Calculate linear optical response (LR-TD-DFTB).
        For details, see Niehaus et.al. Phys. Rev. B 63, 085108 (2001)
                        
        parameters:
        ===========
        calc:         calculator object
        energy_cut:   max energy (in eV) for particle-hole excitations
                      Used to select all particle-hole excitations in the 
                      construction of the matrix, that have excitation energy 
                      less than this value. This implies, that if we are 
                      interested in optical response up to some energy, 
                      energy_cut should be slightly larger.
        timing:       output timing summary after calculation
        out:          output object (file name or object)
        """
        self.calc=calc
        self.st=calc.st
        self.el=calc.el
        self.es=calc.st.es
        self.energy_cut=energy_cut/Hartree
        #self.noc=self.st.get_hoc()+1 #number of occupied states (not index)
        # do not use HOC
        self.nel=self.el.get_number_of_electrons()
        self.norb=self.el.get_nr_orbitals()
        self.e=self.st.get_eigenvalues()[0,:]
        self.f=self.st.get_occupations()[0,:]
        if nu.any( abs(self.st.wf.flatten().imag)>1E-10 ):
            raise ValueError('Wave functions should not be complex.')
        self.wf = self.st.wf[0].real
        self.S = self.st.S[0].real
        self.N=len(self.el)
        if self.calc.get('SCC')==False:
            raise AssertionError('SCC should be True. (Otherwise, just plot DOS)')
        atoms=calc.get_atoms()
        if atoms.pbc.any():
            raise AssertionError('No linear response for extended, periodic systems!')
        
        #if abs(nu.mod(self.nel,2))>1E-2:
            #raise RuntimeError('Linear response only for closed shell systems! (even number of electrons)')
        #if abs(self.nel-2*self.noc)>1E-2:
            #print 'Number of electrons:',self.nel
            #print '2*Number of occupied states:',2*self.noc
            #raise RuntimeError('Number of electrons!=2*number of occupied orbitals. Decrease electronic temperature?')
        if txt is None:
            self.txt=sys.stdout
        else:
            self.txt=open(txt,'a')
        self.timer=Timer('Linear Response',txt=self.txt)   
        self.timing=timing  
        self.done=False       
        self.allowed_cut = 1E-2 #if osc.strength is smaller, transition is not allowed
        self._initialize()
        
    def _initialize(self):
        """
        Perform some initialization calculations.
        """
        self.timer.start('init')
        self.Swf = nu.dot(self.S,self.wf.transpose()) #try to avoind the transpose
        self.timer.stop('init')
            
        
    def get_linear_response(self):
        """ Get linear response spectrum in eV. """
        return self.omega*Hartree, self.F      
        
        
    def mulliken_transfer(self,k,l):
        """ Return Mulliken transfer charges between states k and l. """
        q=[]
        for i,o1,no in self.el.get_property_lists(['i','o1','no']):
            #qi=sum( [self.wf[a,k]*self.Swf[a,l]+self.wf[a,l]*self.Swf[a,k] for a in range(o1,o1+no)] )
            qi=sum( [self.wf[k,a]*self.Swf[a,l]+self.wf[l,a]*self.Swf[a,k] for a in range(o1,o1+no)] )
            q.append(qi/2)
        return nu.array(q)
        
    def run(self):
        """ Run the calculation. """
        if self.done==True:
            raise AssertionError('Run LR calculation only once.')
        
        print>>self.txt, '\nLR for %s (charge %.2f). ' %(self.el.get_name(),self.calc.get_charge()),
        
        #
        # select electron-hole excitations (i occupied, j not occupied)
        # de = excitation energy ej-ei (ej>ei)
        # df = occupation difference fi-fj (ej>ei so that fi>fj)
        #
        de=[] 
        df=[] 
        particle_holes=[]
        gammafct = self.es.get_gamma()
        self.timer.start('setup ph pairs')
        for i in range(self.norb):
            for j in range(i+1,self.norb):
                energy=self.e[j]-self.e[i]
                occup=(self.f[i]-self.f[j])/2 #normalize the double occupations (...is this rigorously right?)
                if energy<self.energy_cut and occup>1E-6:
                    assert energy>0 and occup>0
                    particle_holes.append([i,j])
                    de.append(energy)
                    df.append(occup)
        self.timer.stop('setup ph pairs')
        de=nu.array(de)
        df=nu.array(df)
        
        #
        # setup the matrix (gamma-approximation) and diagonalize
        #
        self.timer.start('setup matrix')
        dim=len(de)
        print>>self.txt, 'Dimension %i. ' %dim,
        if not 0<dim<100000:
            raise RuntimeError('Coupling matrix too large or small (%i)' %dim)
        r=self.el.get_positions()            
        transfer_q=nu.array([self.mulliken_transfer(ph[0],ph[1]) for ph in particle_holes])
        rv=nu.array([dot(tq,r) for tq in transfer_q])
        
        gamma = self.es.get_gamma().copy()        
        gamma_tq=nu.zeros((dim,self.N))
        for k in range(dim):
            gamma_tq[k,:]=dot(gamma,transfer_q[k,:])            
            
        matrix=nu.zeros((dim,dim))            
        for k1,ph1 in enumerate(particle_holes):
            matrix[k1,k1]=de[k1]**2  
            for k2,ph2 in enumerate(particle_holes):
                coupling=dot(transfer_q[k1,:],gamma_tq[k2,:])
                matrix[k1,k2]+=2*sqrt(df[k1]*de[k1]*de[k2]*df[k2])*coupling
        self.timer.stop('setup matrix')                            
                                                                        
        print>>self.txt, 'coupling matrix constructed. ',
        self.txt.flush()         
        self.timer.start('diagonalize')                                                         
        omega2,eigv=eigh(matrix)
        self.timer.stop('diagonalize')
        print>>self.txt, 'Matrix diagonalized.',
        self.txt.flush()
        assert nu.all(omega2>1E-16)           
        omega=sqrt(omega2)
        
        # calculate oscillator strengths
        F=[]
        collectivity=[]
        self.timer.start('oscillator strengths')
        for ex in range(dim):
            v=[]    
            for i in range(3):
                v.append( sum( rv[:,i]*sqrt(df[:]*de[:])*eigv[:,ex])/sqrt(omega[ex])*2 )
            F.append( omega[ex]*dot(v,v)*2.0/3 )
            collectivity.append( 1/sum(eigv[:,ex]**4) )
        self.omega=omega
        self.F=F   
        self.eigv=eigv
        self.collectivity=collectivity
        self.dim=dim
        self.particle_holes=particle_holes
        self.timer.stop('oscillator strengths')
        if self.timing:
            self.timer.summary()
        self.done=True         
        self.emax=max(omega)
        self.particle_holes=particle_holes   
        
            
    def info(self):
        """ Some info about excitations (energy, main p-h excitations,...) """
        print '\n#e(eV), f, collectivity, transitions ...'
        for ex in range(self.dim):
            if self.F[ex]<self.allowed_cut:
                continue
            print '%.5f %.5f %8.1f' %(self.omega[ex]*Hartree,self.F[ex],self.collectivity[ex]), 
            order=nu.argsort(abs(self.eigv[:,ex]))[::-1]
            for ph in order[:4]:
                i,j=self.particle_holes[ph]
                print '%3i-%-3i:%-10.3f' %(i,j,self.eigv[ph,ex]**2),
            print         
            
    def get_excitation(self,i,allowed=True):
        """ Return energy (eV) and oscillation strength for i'th allowed excitation index.
        
        i=0 means first excitation
        """
        if allowed==False:
            return self.omega[i]*Hartree, self.F[i]
        else:
            p=-1
            for k in range(self.dim):
                if self.F[k]>=self.allowed_cut: p+=1
                if p==i:
                    return self.omega[k]*Hartree, self.F[k]                    
                
            
    def write_spectrum(self,filename=None):
        """ Write the linear response spectrum into file. """
        if filename==None:
            filename='linear_spectrum.out'
        o=open(filename,'w')
        print>>o, '#e(eV), f'
        for ex in range(self.dim):
            print>>o, '%10.5f %10.5f %10.5f' %(self.omega[ex]*Hartree,self.F[ex],self.collectivity[ex])
        o.close()  
                  
    def read_spectrum(self,filename):
        """ Read the linear response from given file.
        
        Format: energy & oscillator strength.
        """
        o=open(filename,'r')
        data=mix.read(filename)
        self.omega, self.F, self.collectivity=data[:,0], data[:,1], data[:,2]
                          
                  
    def plot_spectrum(self,filename,width=0.2,xlim=None):
        """ Make pretty plot of the linear response. 
        
        Parameters:
        ===========
        filename: output file name (&format, supported by matplotlib)
        width:    width of Lorenzian broadening 
        xlim:     energy range for plotting tuple (emin,emax)
        """
        import pylab as pl
        if not self.done:
            self.run()
        
        e,f=mix.broaden(self.omega*Hartree,self.F,width=width,N=1000,extend=True) 
        f=f/max(abs(f))
        
        pl.plot(e,f,lw=2)
        xs, ys = pl.poly_between(e, 0, f)
        pl.fill(xs,ys,fc='b',ec='b',alpha=0.5)
        pl.ylim(0,1.2)
        
        if xlim==None:
            pl.xlim(0,self.emax*Hartree*1.2)
        else:
            pl.xlim(xlim)            
        pl.xlabel('energy (eV)')
        pl.ylabel('linear optical response')
        pl.title('Optical response')
        pl.savefig(filename)
        #pl.show()
        pl.close()                  
                
                
                
                
                    
        
        
        
        
        


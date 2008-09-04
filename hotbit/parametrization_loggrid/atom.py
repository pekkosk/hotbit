import numpy as nu
import pylab as pl
from scipy.integrate import odeint
from box.data import atom_occupations
from box.data import atom_valence
from box.data import data
from copy import copy
from box.interpolation import Function 
import sys
from box.misc import AndersonMixer

class KSAllElectron:
    def __init__(self,symbol,occu={},valence=[],confinement=None,etol=1E-8,grid=None,xc='PW92',mixer=None,convergence=1E-7,out=None,itmax=200,verbose=False):
        """ 
        KS_Atom('H',occu={'1s':1},confinement={'type':'frauenheim','r0':1.85},valence=['1s'] 
        
        occu: e.g. {'2s':2,'2p':2}. Overrides (for given orbitals) default occupations and valence from box.data.
        valence: list of valence orbitals, e.g. ['2s','2p']. Overrides default occupations from box.data.
        etol: sp energy tolerance for eigensolver
        grid: RadialGrid class object. Default: ...
        convergence: if integrated |density-density_old|<convergence, density is converged
        itmax: maximum number of iterations for self-consistency.
        mixer: density mixer object (default AndersonMixer with mix=0.3 and memory=3)
        """
        
        self.symbol=symbol
        self.data=copy(data[self.symbol])
        self.Z=self.data['Z']
        self.valence=copy(atom_valence[self.symbol])
        if valence!=[]:
            self.valence=valence
        self.occu=atom_occupations[self.symbol].copy()
        nel_neutral=sum([self.occu[key] for key in self.occu])
        self.occu.update(occu)
        self.nel=sum([self.occu[key] for key in self.occu])
        self.charge=nel_neutral-self.nel
        self.confinement=confinement
        if self.confinement==None:
            self.confinement=ConfinementPotential('none')
        self.convergence=convergence
        self.itmax=itmax      
        if mixer==None:    
            self.mixer=AndersonMixer(0.2,3,1E-50)
        else:
            self.mixer=mixer
        self.verbose=verbose
        self.maxl=9
        self.maxn=9
        self.etol=etol
        self.grid=grid
        self.plotr={}
        self.unl={}
        self.Rnl={}
        self.unl_int={}
        self.Rnl_int={}
        self.xc=xc
        self.total_energy=0.0
        if self.xc=='PW92':
            self.xcf=XC_PW92()
        else:
            raise NotImplementedError('Not implemented: %s' %xc)
        if self.grid==None:
            self.grid=RadialGrid('logarithmic',rmin=1E-6,rmax=120,N=1000)
            self.rgrid=self.grid.get_grid()
        self.nr=self.grid.get_N()
        self.set_output(out)  
        
    def __getstate__(self):
        """ Return dictionary of all pickable items. """
        d=self.__dict__.copy()
        for key in self.__dict__:
            if callable(d[key]):
                d.pop(key)
        d.pop('out')
        return d
        
    def set_output(self,out):
        """ Set output channel and give greetings. """
        if out==None:
            self.out=sys.stdout
        else:
            self.out=open(out,'a')
        print>>self.out, '-------------------------------------------'
        print>>self.out, 'Kohn-Sham all-electron calculation for %2s ' %self.symbol
        print>>self.out, '-------------------------------------------'
        print>>self.out, self.get_comment()
        
    def calculate_energies(self,echo=False):
        """ Calculate energy contributions. """
        self.bs_energy=0.0
        for n,l,nl in self.list_states():
            self.bs_energy+=self.occu[nl]*self.enl[nl]
        self.Hartree_energy=self.grid.integrate(self.Hartree*self.dens,use_dV=True)/2
        self.vxc_energy=self.grid.integrate(self.vxc*self.dens,use_dV=True)
        self.exc_energy=self.grid.integrate(self.exc*self.dens,use_dV=True)
        self.confinement_energy=self.grid.integrate(self.conf*self.dens,use_dV=True)
        self.total_energy=self.bs_energy-self.Hartree_energy-self.vxc_energy+self.exc_energy
        if echo:
            print>>self.out, 'single-particle energies'
            print>>self.out, '------------------------'
            for n,l,nl in self.list_states():
                print>>self.out, 'nl, energy %.15f' %self.enl[nl]
            print>>self.out, '\n\n'
            print>>self.out, 'energies:'
            print>>self.out, '---------'
            print>>self.out, 'sum of eigenvalues:     %.15f' %self.bs_energy
            print>>self.out, 'Hartree energy:         %.15f' %self.Hartree_energy
            print>>self.out, 'vxc correction:         %.15f' %self.vxc_energy
            print>>self.out, 'exchange + corr energy: %.15f' %self.exc_energy
            print>>self.out, '----------------------------'
            print>>self.out, 'total energy:           %.15f' %self.total_energy     
        
    def calculate_density(self):
        """ Calculate the radial electron density.; sum_nl |Rnl(r)|**2/(4*pi) """
        dens=nu.zeros_like(self.rgrid)
        for n,l,nl in self.list_states():
            dens+=self.occu[nl]*self.unl[nl]**2/(4*nu.pi*self.rgrid**2)
        self.dens=dens
        self.nel_error=self.grid.integrate(dens,use_dV=True)-self.nel
      
    def calc_Hartree(self):
        """ Calculate Hartree potential. """
        dV=self.grid.get_dvolumes()
        r0=self.grid.get_r0grid()
        nr=self.grid.get_N()
        n0=0.5*(self.dens[1:nr]+self.dens[0:nr-1])
        hi=(dV*n0/r0).sum()
        lo=0.0
        Hartree=nu.zeros(nr)
        Hartree[0]=hi
        for i in range(nr-1):
            lo+=dV[i]*n0[i]
            hi-=dV[i]*n0[i]/r0[i]
            Hartree[i+1]=lo/r0[i]+hi
        self.Hartree=Hartree
    
    def V_nuclear(self,r):
        return -self.Z/r
    
    def calc_veff(self):
        """ Calculate effective potential. """
        self.veff=nu.zeros_like(self.rgrid)
        self.vxc=nu.zeros_like(self.rgrid)
        self.exc=nu.zeros_like(self.rgrid)
        self.conf=nu.zeros_like(self.rgrid)
        for i,r in enumerate(self.rgrid):
            self.vxc[i]=self.xcf.vxc(self.dens[i])
            self.exc[i]=self.xcf.exc(self.dens[i])
            self.conf[i]=self.confinement(r)
            self.veff[i]=self.V_nuclear(r) + self.Hartree[i] + self.vxc[i] + self.conf[i]
            
        self.V_effective=Function('spline',self.rgrid,self.veff)
        #self.V_effective=FastSpline(self.rgrid,self.veff,'power',p=3,xmin=1E-6,xmax=100.0)
        #self.V_effective=Function('fastspline',self.rgrid,self.veff,'power',p=3,xmin=1E-6,xmax=100.0)
        #self.V_effective.plot()
        
    
    def solve_ground_state(self):
        """
        Solve the self-consistent potential. 
        """
        self.enl={}
        self.d_enl={}
        for n,l,nl in self.list_states():
            self.enl[nl]=0.0
            self.d_enl[nl]=0.0
        self.rgrid=self.grid.get_grid() 
        self.nr=len(self.rgrid)
        
        N=self.grid.get_N()
        self.dens=nu.zeros((N,))
        self.Hartree=nu.zeros((N,))
        for it in range(self.itmax):
            self.calc_veff()
            self.solve_eigenstates2(it)
            
            dens0=self.dens.copy()
            total_energy0=self.total_energy
            self.calculate_density()
            diff=self.grid.integrate(nu.abs(self.dens-dens0),use_dV=True)
            dnel=abs(self.grid.integrate(self.dens,use_dV=True)-self.nel)
            done,self.dens=self.mixer(dens0,self.dens)
            if diff<self.convergence: # and dnel<self.convergence:
                break
            self.calc_Hartree()
            self.calculate_energies()
            ediff=abs(total_energy0-self.total_energy)
            print>>self.out, 'iteration %3i, energy %9.5f, de %9.6f (%.3g,%.3g<%.3g)' \
                            %(it,self.total_energy,ediff,diff,dnel,self.convergence)
            if it==self.itmax-1:
                raise RuntimeError('Density not converged in %i iterations' %(it+1))
        
        #print self.grid.integrate(self.unl['1s']**2)
        #print>>self.out, 'converged in %i iterations.' %(it+1)
        print>>self.out, '%9.4f electrons, %9.4f should be' %(self.grid.integrate(self.dens,use_dV=True),self.nel)
        self.calculate_energies(echo=True)
        for n,l,nl in self.list_states():
            self.Rnl_int[nl]=Function('spline',self.rgrid,self.Rnl[nl])
            self.unl_int[nl]=Function('spline',self.rgrid,self.unl[nl])
                  
    def __function(self,y,rho,l,epsilon):
        """ The integration function needed for odeint. """
        r=nu.exp(rho)
        return nu.array( [y[1],y[1]+2*(l*(l+1)+r**2*(self.V_effective(r)-epsilon))*y[0]] )
        #return nu.array( [y[1],2*(0.5*l*(l+1)+r**2*(self.V_effective(r)-epsilon)-0.25)*y[0]] )
    
    def solve_eigenstates2(self,iteration,itmax=100):
        """ 
        Solve the eigenstates for given effective potential.
        
        """
        M=100
        grid=nu.linspace(0,6,M)
        for n,l,nl in self.list_states():
            nodes_nl=n-l-1
            eps=self.guess_epsilon(init=True,iteration=iteration,nl=nl)   
            u0=nu.exp(grid[0]*(l+1))
            dudr0=(l+1)*nu.exp(grid[0]*(l+1))
            for it in range(itmax):
                u=odeint(self.__function,[u0,dudr0],grid,args=(l,eps))
                nodes=sum( (u[0:M-1,0]*u[1:M,0])<0 )
                print eps,nodes
                eps,done=self.guess_epsilon(eps,nodes,u[-1,0],nodes_nl,nl)
                if done:
                    break
            
            
            
            rg=nu.exp(grid)
            
            self.unl[nl]=self.normalize(nl,u[:,0])
            self.Rnl[nl]=self.unl[nl]/self.rgrid
            self.d_enl[nl]=abs(eps-self.enl[nl])
            self.enl[nl]=eps
            if self.verbose:
                print>>self.out, '-- state %s, %i eigensolver iterations, e=%9.5f, de=%9.5f' %(nl,it,self.enl[nl],self.d_enl[nl])
                                  
              
    def _function(self,y,r,l,epsilon):
        """ The integration function needed for odeint. """
        return nu.array( [y[1],2*(self.V_effective(r)+l*(l+1)/(2*r**2)-epsilon)*y[0]] )
    
    def solve_eigenstates(self,iteration,itmax=100):
        """ 
        Solve the eigenstates for given effective potential.
        
        For small r u(r)~r**(l+1). Then u(r0)=C*r0**(l+1) and u'(r0)=C*(l+1)*r0**l; C comes from normalization.
        """
        for n,l,nl in self.list_states():
            nodes_nl=n-l-1
            eps=self.guess_epsilon(init=True,iteration=iteration,nl=nl)   
            u0=self.rgrid[0]**(l+1)
            dudr0=(l+1)*self.rgrid[0]**l                
            for it in range(itmax):
                try:
                    u=odeint(self._function,[u0,dudr0],self.rgrid,args=(l,eps))
                except:
                    raise RuntimeError
                nodes=sum( (u[0:self.nr-1,0]*u[1:self.nr,0])<0 )
                eps,done=self.guess_epsilon(eps,nodes,u[-1,0],nodes_nl,nl)
                if done:
                    break
            
            self.unl[nl]=self.normalize(nl,u[:,0])
            self.Rnl[nl]=self.unl[nl]/self.rgrid
            self.d_enl[nl]=abs(eps-self.enl[nl])
            self.enl[nl]=eps
            if self.verbose:
                print>>self.out, '-- state %s, %i eigensolver iterations, e=%9.5f, de=%9.5f' %(nl,it,self.enl[nl],self.d_enl[nl])
        
    def guess_epsilon(self,epsilon=None,nodes=None,yn=None,nodes_nl=None,nl=None,init=False,iteration=None):
        """ Guess next single-particle energy. """
        if init:
            if iteration<3:
                (n,l)=orbit_transform(nl,string=False)
                self.high=[self.Z**2/n**2,100,None]
                self.low=[-self.high[0],-1,None]
                eps=-0.5*self.Z**2/n**2
            else:
                eps=self.enl[nl]
                self.high=[self.enl[nl]+2*self.d_enl[nl],100,None]
                self.low=[self.enl[nl]-2*self.d_enl[nl],0,None]
                
            return eps
        elif epsilon!=None:    
            if nodes>nodes_nl:
                self.high=[epsilon,nodes,yn]
            elif nodes<=nodes_nl:
                self.low=[epsilon,nodes,yn]
            if self.high[1]<=self.low[1]:
                raise RuntimeError('Bracketing eigenvalue failed. Use larger tolerances.')
            #if self.high[2]!=None and self.low[2]!=None and self.high[1]-self.high[1]==1:
                #eps=self.low[0]-self.low[2]*(self.high[0]-self.low[0])/(self.high[2]-self.low[2])
            #else:
            eps=0.5*(self.low[0]+self.high[0])
            done=False
            if self.high[0]-self.low[0]<self.etol:
                assert self.high[1]-self.low[1]==1
                done=True
            return eps,done
        
    def normalize(self,nl,u):
        """ 
        Normalize the wave function u; int_r |u(r)|**2 dr=1. Normalize also
        in the sense of cutting the crazy tails with large r.
        """
        mx=1E900
        for i in range(self.nr)[::-1]:
            if abs(u[i])<mx:
                mx=abs(u[i])
            else:
                u[i:]=0.0
                break
            
        # search convenient range for plotting
        Rnl=u/self.rgrid            
        maxval=abs(Rnl).max()
        for i in range(self.nr)[::-1]:            
            if abs(Rnl[i])>1E-3*maxval:
                self.plotr[nl]=i
                break
            
        norm=self.grid.integrate(u**2)
        u=u*nu.sign(u[5]) #sign convention: first antinode positive
        assert norm>1E-6 and norm<1E500
        return u/nu.sqrt(norm)
        
    def plot_Rnl(self,screen=False):
        """ Plot radial wave functions with matplotlib. """
        i=1
        for n,l,nl in self.list_states():
            pl.subplot(3,3,i)
            ri=self.plotr[nl]
            pl.plot(self.rgrid[:ri],self.Rnl[nl][:ri])
            pl.xlabel('r (Bohr)')
            pl.ylabel('Rnl (%s)' %nl)
            i+=1
        if screen:
            pl.show()
        else:
            pl.savefig('%s_KSatom.png' %self.symbol)
        
    def list_states(self):
        """ List all potential states {(n,l,'nl')}. """
        states=[]
        for l in range(self.maxl+1):
            for n in range(1,self.maxn+1):  
                nl=orbit_transform((n,l),string=True) 
                if nl in self.occu:
                    states.append((n,l,nl))
        return states
                    
    def get_eigenvalue(self,nl):
        """ get_eigenvalue('2p') or get_eigenvalue((2,1)) """
        nls=orbit_transform(nl,string=True)
        return self.enl[nls]
        
    def Rnl(self,r,nl):
        """ Rnl(r,'2p') or Rnl(r,(2,1))"""
        nls=orbit_transform(nl,string=True)
        return self.Rnl_int[nls](r)
        
    def unl(self,r,nl):
        """ unl(r,'2p')=Rnl(r,'2p')/r or unl(r,(2,1))..."""
        nls=orbit_transform(nl,string=True)
        return self.unl_int[nls](r)
        
    def get_valence(self):
        """ Get list of valence orbitals, e.g. ['2s','2p'] """
        return self.valence
    
    def get_symbol(self):
        """ Return atom's chemical symbol. """
        return self.symbol
        
    def get_comment(self):
        """ One-line comment, e.g. 'H, charge=0, frauenheim, r0=4' """
        comment='%s xc=%s charge=%.1f conf:%s' %(self.symbol,self.xc,float(self.charge),self.confinement.get_comment())
        return comment
    
    def get_valence_energies(self):
        """ Return list of valence energies, e.g. ['2s','2p'] --> [-39.2134,-36.9412] """
        return [(nl,self.enl[nl]) for nl in self.valence]
    
    def get_valence_l(self):
        pass
    
    def write_functions(self,file,only_valence=True):
        """ Write functions (unl,v_effective,...) into file (only valence functions by default). """
        if only_valence:
            orbitals=self.valence
        else:
            orbitals=[nl for n,l,nl in self.list_states()]
        o=open(file,'a')
        for nl in orbitals:
            print>>o, '\n\nu_%s=' %nl
            for r,u in zip(self.rgrid,self.unl[nl]):
                print>>o, r,u
            
        print>>o,'\n\nv_effective='
        for r,ve in zip(self.rgrid,self.veff):
                print>>o, r,ve        
        print>>o,'\n\nconfinement='
        for r,vc in zip(self.rgrid,self.conf):
                print>>o, r,vc
        print>>o,'\n\n'
        
    
class RadialGrid:
    def __init__(self,mode,**kwargs):
        """ 
        mode
        ----
        'linear' r(i)=(i-1)/(N-1)*rmax (kwargs: 'rmax','N')
        'quadratic' r(i)=rmin + (rmax-rmin)/(N-1)**2*i**2 (i=0,...,N-1)
        'exponential' r(i)=rmin*(rmax/rmin)**(i/N) (kwargs: 'rmin','rmax','N')
        
        
        rmin                                                        rmax
        r[0]     r[1]      r[2]            ...                     r[N-1] grid
        I----'----I----'----I----'----I----'----I----'----I----'----I
           r0[0]     r0[1]     r0[2]       ...              r0[N-2]       r0grid
           dV[0]     dV[1]     dV[2]       ...              dV[N-2]       dV
           
           dV[i] is volume element of shell between r[i] and r[i+1]
        """
        self.mode=mode
        self.kwargs=kwargs
        self.make_grid()

    def make_grid(self):
        rmin=self.kwargs['rmin']
        rmax=self.kwargs['rmax']
        N=self.kwargs['N']
        self.rmin,self.rmax,self.N=rmin,rmax,N
        if self.mode=='linear':
            self.grid=nu.linspace(rmin,rmax,N)
        elif self.mode=='exponential':
            self.grid=rmin*(rmax/rmin)**(nu.linspace(1,N,N)/N)        
        elif self.mode=='quadratic':
            self.grid=rmin+(rmax-rmin)/(N-1.0)**2*nu.linspace(0,N-1,N)**2
        elif self.mode=='logarithmic':
            self.grid=rmin*(rmax/rmin)**(nu.arange(N)/(N-1.0))
        elif self.mode=='qubic':
            self.grid=rmin+(rmax-rmin)/(N-1.0)**3*nu.linspace(0,N-1,N)**3
        self.dr=self.grid[1:N]-self.grid[0:N-1]
        self.r0=self.grid[0:N-1]+self.dr/2
        # first dV is sphere, others are shells
        self.dV=nu.zeros_like(self.r0)
        self.dV[0]=4.0/3*nu.pi*self.grid[0]**3 
        self.dV[1:N]=4.0/3*nu.pi*(self.grid[2:N]**3-self.grid[1:N-1]**3)
        return self.grid

    def get_grid(self):
        """ Return the whole radial grid. """
        return self.grid

    def get_N(self):
        """ Return the number of grid points. """
        return self.N
        
    def get_drgrid(self):
        """ Return the grid spacings (array of length N-1). """
        return self.dr
        
    def get_r0grid(self):
        """ Return the mid-points between grid spacings (array of length N-1). """
        return self.r0
    
    def get_dvolumes(self):
        """ Return dV(r)'s=4*pi*r**2*dr. """
        return self.dV
        
    def plot(self,screen=True):
        rgrid=self.get_grid()
        pl.scatter(range(len(rgrid)),rgrid)
        if screen: pl.show()       
         
    def integrate(self,f,use_dV=False):
        """ 
        Integrate function f (given with N grid points).
        int_rmin^rmax f*dr (use_dv=False) or int_rmin^rmax*f dV (use_dV=True)
        """
        if use_dV:
            return ((f[0:self.N-1]+f[1:self.N])*self.dV).sum()*0.5
        else:
            return ((f[0:self.N-1]+f[1:self.N])*self.dr).sum()*0.5
        
        
        
        
class ConfinementPotential:
    def __init__(self,mode,**kwargs):
        self.mode=mode
        if mode=='none':
            self.f=self.none #lambda r:0.0
            self.comment='none'
        elif mode=='frauenheim':
            self.r0=kwargs['r0']
            self.f=self.frauenheim #lambda r:(r/self.r0)**2
            self.comment='frauenheim r0=%.3f' %self.r0
        else:
            raise NotImplementedError('implement new confinements')
        
    def none(self,r):
        return 0.0
        
    def frauenheim(self,r):
        return (r/self.r0)**2
        
    def __call__(self,r):
        return self.f(r)
        
    def get_comment(self):
        return self.comment
        






class XC_PW92:
    def __init__(self):
        self.small=1E-90
        self.a1 = 0.21370
        self.c0 = 0.031091
        self.c1 = 0.046644
        self.b1 = 1.0/2.0/self.c0*nu.exp(-self.c1/2.0/self.c0)
        self.b2 = 2*self.c0*self.b1**2
        self.b3 = 1.6382
        self.b4 = 0.49294
    
    def exc(self,n):
        """ Exchange-correlation with electron density n. """
        if n<self.small:
            return 0.0
        else:
            return self.e_x(n)+self.e_corr(n)
        
    def e_x(self,n):
        """ Exchange. """
        return -3.0/4*(3*n/nu.pi)**(1.0/3)

    def e_corr(self,n):
        """ Correlation energy. """
        rs = (3.0/(4*nu.pi*n))**(1.0/3)
        ec = -2*self.c0*(1+self.a1*rs)*nu.log(1+1/(2*self.c0*(self.b1*nu.sqrt(rs)+self.b2*rs+self.b3*rs**(3.0/2)+self.b4*rs**2)))
        return ec
    
    def vxc(self,n):
        """ Exchange-correlation potential (functional derivative of exc). """
        eps=1E-9*n       
        if n<self.small:
            return 0.0
        else: 
            return ((n+eps)*self.exc(n+eps)-(n-eps)*self.exc(n-eps)) / (2*eps)


    
    
angular_momenta=['s','p','d','f','g','h','i','j','k','l']
def orbit_transform(nl,string):
    """ Transform orbitals into strings<->tuples, e.g. (2,1)<->'2p'. """
    if string==True and type(nl)==type(''):
        return nl #'2p'->'2p'
    elif string==True:
        return '%i%s' %(nl[0],angular_momenta[nl[1]]) #(2,1)->'2p'
    elif string==False and type(nl)==type((2,1)):
        return nl      #(2,1)->(2,1)
    elif string==False:
        l=angular_momenta.index(nl[1])
        n=int(nl[0])
        return (n,l)  # '2p'->(2,1)
        
            
    


#occupations={'H':{'1s':1}}
#valence={'H':['1s']}

  
#def get_density_FWHM
#def calculate_ionization_energy
#def calculate_electron_affinity
#def calculate_hubbard_U


 
        
    
    
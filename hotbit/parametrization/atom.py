import numpy as nu
import pylab as pl
from scipy.integrate import odeint
from box.data import data
from copy import copy
from box.interpolation import Function, SplineFunction
import sys
import os
from box.timing import Timer
from time import asctime
import math
array=nu.array
sqrt=math.sqrt
pi=math.pi
log=math.log

class KSAllElectron:
    def __init__(self,symbol,
                      configuration={},
                      valence=[],
                      confinement=None,
                      xc='PW92',
                      convergence={'density':1E-7,'energies':1E-7},        
                      ScR=False,
                      nodegpts=500,
                      mix=0.2,
                      itmax=200,  
                      timing=False,
                      verbose=False,
                      txt=None):       
        """ 
        Make Kohn-Sham all-electron calculation for given atom.
        
        Examples:
        ---------
        atom=KSAllElectron('C')
        atom=KSAllElectron('C',confinement={'mode':'frauenheim','r0':1.234})
        atom.run()
        
        Parameters:
        -----------
        symbol:         chemical symbol
        configuration:  e.g. {'2s':2,'2p':2}. Overrides (for orbitals given in dict) default 
                        configuration from box.data.
        valence:        valence orbitals, e.g. ['2s','2p']. Overrides default valence from box.data.
        confinement:    additional confining potential (see ConfinementPotential class)
        etol:           sp energy tolerance for eigensolver (Hartree)
        convergence:    convergence criterion dictionary
                        * density: max change for integrated |n_old-n_new| 
                        * energies: max change in single-particle energy (Hartree)
        ScR:            Use scalar relativistic corrections
        nodegpts:       total number of grid points is nodegpts times the max number 
                        of antinodes for all orbitals
        mix:            effective potential mixing constant
        itmax:          maximum number of iterations for self-consistency.
        timing:         output of timing summary 
        verbose:        increase verbosity during iterations
        txt:            output file name 
        """
        self.symbol=symbol
        self.valence=valence
        self.confinement=confinement
        self.xc=xc
        self.convergence=convergence
        self.ScR = ScR
        self.set_output(txt)
        self.itmax=itmax
        self.verbose=verbose
        self.nodegpts=nodegpts
        self.mix=mix
        self.timing=timing       
        self.timer=Timer('KSAllElectron',txt=self.txt,enabled=self.timing)
        self.timer.start('init')
        
        # element data
        self.data=copy( data[self.symbol] )
        self.Z=self.data['Z']
        if self.valence == []:
            self.valence = copy( data[self.symbol]['valence_orbitals'] )
        
        # ... more specific
        self.occu = copy( data[self.symbol]['configuration'] )
        nel_neutral = self.Z
        assert sum([self.occu[key] for key in self.occu]) == nel_neutral
        self.occu.update( configuration )
        self.nel=sum([self.occu[key] for key in self.occu])
        self.charge=nel_neutral-self.nel
        if self.confinement==None:
            self.confinement_potential=ConfinementPotential('none')
        else:
            self.confinement_potential=ConfinementPotential(**self.confinement)            
        self.conf=None      
        self.nucl=None      
        self.exc=None
        if self.xc=='PW92':
            self.xcf=XC_PW92()
        else:
            raise NotImplementedError('Not implemented xc functional: %s' %xc)
        
        # technical stuff
        self.maxl=9
        self.maxn=9
        self.plotr={}
        self.unlg={}
        self.Rnlg={}
        self.unl_fct={}
        self.Rnl_fct={}
        self.veff_fct=None
        self.total_energy=0.0
        
        maxnodes=max( [n-l-1 for n,l,nl in self.list_states()] )        
        self.rmin, self.rmax, self.N=( 1E-2/self.Z, 100.0, (maxnodes+1)*self.nodegpts )
        if self.ScR:
            print >> self.txt, 'Using scalar relativistic corrections.'
        print>>self.txt, 'max %i nodes, %i grid points' %(maxnodes,self.N)
        self.xgrid=nu.linspace(0,nu.log(self.rmax/self.rmin),self.N)
        self.rgrid=self.rmin*nu.exp(self.xgrid)
        self.grid=RadialGrid(self.rgrid)
        self.timer.stop('init')
        print>>self.txt, self.get_comment()
        self.solved=False
        
    def __getstate__(self):
        """ Return dictionary of all pickable items. """
        d=self.__dict__.copy()
        for key in self.__dict__:
            if callable(d[key]):
                d.pop(key)
        d.pop('out')
        return d
        
    def set_output(self,txt):
        """ Set output channel and give greetings. """
        if txt == '-':
            self.txt = open(os.devnull,'w')
        elif txt==None:
            self.txt=sys.stdout
        else:
            self.txt=open(txt,'a')
        print>>self.txt, '*******************************************'
        print>>self.txt, 'Kohn-Sham all-electron calculation for %2s ' %self.symbol
        print>>self.txt, '*******************************************'
    
    
    def calculate_energies(self,echo=False):
        """ 
        Calculate energy contributions. 
        """
        self.timer.start('energies')
        self.bs_energy=0.0
        for n,l,nl in self.list_states():
            self.bs_energy+=self.occu[nl]*self.enl[nl]
            
        self.exc=array([self.xcf.exc(self.dens[i]) for i in xrange(self.N)])
        self.Hartree_energy=self.grid.integrate(self.Hartree*self.dens,use_dV=True)/2
        self.vxc_energy=self.grid.integrate(self.vxc*self.dens,use_dV=True)
        self.exc_energy=self.grid.integrate(self.exc*self.dens,use_dV=True)
        self.confinement_energy=self.grid.integrate(self.conf*self.dens,use_dV=True)
        self.total_energy=self.bs_energy-self.Hartree_energy-self.vxc_energy+self.exc_energy
        if echo:
            print>>self.txt, '\n\nEnergetics:'
            print>>self.txt, '-------------'
            print>>self.txt, '\nsingle-particle energies'
            print>>self.txt, '------------------------'
            for n,l,nl in self.list_states():
                print>>self.txt, '%s, energy %.15f' %(nl,self.enl[nl])
                            
            print>>self.txt, '\nvalence orbital energies'
            print>>self.txt, '--------------------------'
            for nl in data[self.symbol]['valence_orbitals']:
                print>>self.txt, '%s, energy %.15f' %(nl,self.enl[nl])                            
                
            print>>self.txt, '\n'
            print>>self.txt, 'total energies:'
            print>>self.txt, '---------------'
            print>>self.txt, 'sum of eigenvalues:     %.15f' %self.bs_energy
            print>>self.txt, 'Hartree energy:         %.15f' %self.Hartree_energy
            print>>self.txt, 'vxc correction:         %.15f' %self.vxc_energy
            print>>self.txt, 'exchange + corr energy: %.15f' %self.exc_energy
            print>>self.txt, '----------------------------'
            print>>self.txt, 'total energy:           %.15f\n\n' %self.total_energy     
        self.timer.stop('energies')            
        
        
    def calculate_density(self):
        """ Calculate the radial electron density.; sum_nl |Rnl(r)|**2/(4*pi) """
        self.timer.start('density')
        dens=nu.zeros_like(self.rgrid)
        for n,l,nl in self.list_states():
            dens+=self.occu[nl]*self.unlg[nl]**2
                   
        nel=self.grid.integrate(dens)
        if abs(nel-self.nel)>1E-10:
            raise RuntimeError('Integrated density %.3g, number of electrons %.3g' %(nel,self.nel) )
        dens=dens/(4*nu.pi*self.rgrid**2)
        
        self.timer.stop('density')
        return dens           
        
      
    def calculate_Hartree_potential(self):
        """ 
        Calculate Hartree potential. 
        
        Everything is very sensitive to the way this is calculated.
        If you can think of how to improve this, please tell me!
        
        """
        self.timer.start('Hartree')
        dV=self.grid.get_dvolumes()
        r, r0=self.rgrid, self.grid.get_r0grid()
        N=self.N
        n0=0.5*(self.dens[1:]+self.dens[:-1])
        n0*=self.nel/sum(n0*dV)
        
        lo, hi, Hartree=nu.zeros(N), nu.zeros(N), nu.zeros(N)        
        lo[0]=0.0 
        for i in range(1,N):    
            lo[i] = lo[i-1] + dV[i-1]*n0[i-1] 
            
        hi[-1]=0.0
        for i in range(N-2,-1,-1):            
            hi[i] = hi[i+1] + n0[i]*dV[i]/r0[i]
            
        for i in range(N):
            Hartree[i] = lo[i]/r[i] + hi[i]                            
        self.Hartree=Hartree
        self.timer.stop('Hartree')
    
    
    def V_nuclear(self,r):
        return -self.Z/r
    
    
    def calculate_veff(self):
        """ Calculate effective potential. """
        self.timer.start('veff')
        self.vxc=array([self.xcf.vxc(self.dens[i]) for i in xrange(self.N)])
        self.timer.stop('veff')
        return self.nucl + self.Hartree + self.vxc + self.conf        
    
    
    def guess_density(self):
        """ Guess initial density. """
        r2=0.02*self.Z # radius at which density has dropped to half; improve this!
        dens=nu.exp( -self.rgrid/(r2/nu.log(2)) )
        dens=dens/self.grid.integrate(dens,use_dV=True)*self.nel
        pl.plot(self.rgrid,dens)
        return dens
        
    
    def run(self):
        """
        Solve the self-consistent potential. 
        """
        self.timer.start('solve ground state')
        print>>self.txt, '\nStart iteration...'
        self.enl={}
        self.d_enl={}
        for n,l,nl in self.list_states():
            self.enl[nl]=0.0
            self.d_enl[nl]=0.0
         
        N=self.grid.get_N()
        
        # make confinement and nuclear potentials; intitial guess for veff
        self.conf=array([self.confinement_potential(r) for r in self.rgrid])
        self.nucl=array([self.V_nuclear(r) for r in self.rgrid])
        self.veff=self.nucl+self.conf
        
        self.dens=self.guess_density()
        self.calculate_Hartree_potential()
        #self.Hartree=nu.zeros((N,))
        
        for it in range(self.itmax):
            self.veff=self.mix*self.calculate_veff()+(1-self.mix)*self.veff
            if self.ScR:
                veff = SplineFunction(self.rgrid, self.veff)
                self.dveff = array([veff(r, der=1) for r in self.rgrid])
            d_enl_max, itmax=self.solve_eigenstates(it)
            
            dens0=self.dens.copy()
            self.dens=self.calculate_density()
            diff=self.grid.integrate(nu.abs(self.dens-dens0),use_dV=True)
                        
            if diff<self.convergence['density'] and d_enl_max<self.convergence['energies']:
                break
            self.calculate_Hartree_potential()
            if nu.mod(it,10)==0:
                print>>self.txt, 'iter %3i, dn=%.1e>%.1e, max %i sp-iter' %(it,diff,self.convergence['density'],itmax)
            if it==self.itmax-1:
                if self.timing:
                    self.timer.summary()            
                raise RuntimeError('Density not converged in %i iterations' %(it+1))
            self.txt.flush()            
        
        self.calculate_energies(echo=True)        
        print>>self.txt, 'converged in %i iterations' %it
        print>>self.txt, '%9.4f electrons, should be %9.4f' %(self.grid.integrate(self.dens,use_dV=True),self.nel)
        for n,l,nl in self.list_states():
            self.Rnl_fct[nl]=Function('spline',self.rgrid,self.Rnlg[nl])
            self.unl_fct[nl]=Function('spline',self.rgrid,self.unlg[nl])
        self.timer.stop('solve ground state')                    
        self.timer.summary()    
        self.txt.flush()
        self.solved=True


    def solve_eigenstates(self,iteration,itmax=100):
        """ 
        Solve the eigenstates for given effective potential.
        
        u''(r) - 2*(v_eff(r)+l*(l+1)/(2r**2)-e)*u(r)=0 
        ( u''(r) + c0(r)*u(r) = 0 )
        
        r=r0*exp(x) --> (to get equally spaced integration mesh)
        
        u''(x) - u'(x) + c0(x(r))*u(r) = 0
        """
        self.timer.start('eigenstates')

        rgrid=self.rgrid
        xgrid=self.xgrid
        dx=xgrid[1]-xgrid[0]
        N=self.N
        c2=nu.ones(N)
        c1=-nu.ones(N)
        d_enl_max=0.0
        itmax=0

        for n,l,nl in self.list_states():
            nodes_nl=n-l-1
            if iteration==0:
                eps=-1.0*self.Z**2/n**2

            else:
                eps=self.enl[nl]

            if iteration<=3:
                delta=0.5*self.Z**2/n**2  #previous!!!!!!!!!!                
            else:
                delta=self.d_enl[nl]

            direction='none'
            epsmax=self.veff[-1]-l*(l+1)/(2*self.rgrid[-1]**2)
            it=0
            u=nu.zeros(N)
            hist=[]

            while True:
                eps0=eps
                c0, c1, c2 = self.construct_coefficients(l, eps)

                # boundary conditions for integration from analytic behaviour (unscaled)
                # u(r)~r**(l+1)   r->0
                # u(r)~exp( -sqrt(c0(r)) ) (set u[-1]=1 and use expansion to avoid overflows)
                u[0:2]=rgrid[0:2]**(l+1)

                if not(c0[-2]<0 and c0[-1]<0):
                    pl.plot(c0)
                    pl.show()

                assert c0[-2]<0 and c0[-1]<0

                u, nodes, A, ctp=shoot(u,dx,c2,c1,c0,N)
                it+=1
                norm=self.grid.integrate(u**2)
                u=u/sqrt(norm)

                if nodes>nodes_nl:
                    # decrease energy
                    if direction=='up': delta/=2
                    eps-=delta
                    direction='down'
                elif nodes<nodes_nl:
                    # increase energy
                    if direction=='down': delta/=2
                    eps+=delta
                    direction='up'
                elif nodes==nodes_nl:
                    shift=-0.5*A/(rgrid[ctp]*norm)
                    if abs(shift)<1E-8: #convergence
                        break
                    if shift>0:
                        direction='up'
                    elif shift<0:
                        direction='down'
                    eps+=shift
                if eps>epsmax:
                    eps=0.5*(epsmax+eps0)
                hist.append(eps)

                if it>100:
                    print>>self.txt, 'Epsilon history for %s' %nl
                    for h in hist:
                        print h
                    print>>self.txt, 'nl=%s, eps=%f' %(nl,eps)
                    print>>self.txt, 'max epsilon',epsmax
                    raise RuntimeError('Eigensolver out of iterations. Atom not stable?')

            itmax=max(it,itmax)
            self.unlg[nl]=u
            self.Rnlg[nl]=self.unlg[nl]/self.rgrid
            self.d_enl[nl]=abs(eps-self.enl[nl])
            d_enl_max=max(d_enl_max,self.d_enl[nl])
            self.enl[nl]=eps
            if self.verbose:
                print>>self.txt, '-- state %s, %i eigensolver iterations, e=%9.5f, de=%9.5f' %(nl,it,self.enl[nl],self.d_enl[nl])

            assert nodes==nodes_nl
            assert u[1]>0.0
        self.timer.stop('eigenstates')
        return d_enl_max, itmax


    def construct_coefficients(self, l, eps):
        c = 137.036
        c2 = nu.ones(self.N)
        if self.ScR == False:
            c0 = -2*( 0.5*l*(l+1)+self.rgrid**2*(self.veff-eps) )
            c1 = -nu.ones(self.N)
        else:
            # from Paolo Giannozzi: Notes on pseudopotential generation
            ScR_mass = array([1 + 0.5*(eps-V)/c**2 for V in self.veff])
            c0 = -l*(l+1) - 2*ScR_mass*self.rgrid**2*(self.veff-eps) - self.dveff*self.rgrid/(2*ScR_mass*c**2)
            c1 = self.rgrid*self.dveff/(2*ScR_mass*c**2) - 1
        return c0, c1, c2


    def plot_Rnl(self,screen=False,r=False):
        """ Plot radial wave functions with matplotlib. 
        r: plot as a function of r or grid points
        """
        i=1
        rmax=data[self.symbol]['R_cov']/0.529177*2
        ri=self.N #sum(self.rgrid<rmax)
        states=len(self.list_states())
        p=nu.ceil(nu.sqrt(states)) #p**2>=states subplots
        for n,l,nl in self.list_states():
            pl.subplot(p,p,i)
            if r:
                pl.plot(self.rgrid[:ri],self.Rnlg[nl][:ri])    
                pl.xlabel('r (Bohr)')
            else:                
                pl.plot(self.Rnlg[nl])
                pl.xlabel('r (grid point)')
            pl.ylabel('Rnl (%s)' %nl)
            i+=1
        if screen:
            pl.show()
        else:
            pl.savefig('%s_KSatom.png' %self.symbol)
            
            
    def get_wf_range(self,nl,fractional_limit=1E-7):
        """ Return the maximum r for which |R(r)|<fractional_limit*max(|R(r)|) """
        wfmax=max(abs(self.Rnlg[nl]))
        for r,wf in zip(self.rgrid[-1::-1],self.Rnlg[nl][-1::-1]):
            if abs(wf)>fractional_limit*wfmax: 
                return r
        
        
    def list_states(self):
        """ List all potential states {(n,l,'nl')}. """
        states=[]
        for l in range(self.maxl+1):
            for n in range(1,self.maxn+1):  
                nl=orbit_transform((n,l),string=True) 
                if nl in self.occu:
                    states.append((n,l,nl))
        return states
                
                
    def get_energy(self):
        return self.total_energy                
                    
                    
    def get_epsilon(self,nl):
        """ get_eigenvalue('2p') or get_eigenvalue((2,1)) """
        nls=orbit_transform(nl,string=True)
        if not self.solved:
            raise AssertionError('run calculations first.')
        return self.enl[nls]
    
    
    def effective_potential(self,r,der=0):
        """ Return effective potential at r or its derivatives. """
        if self.veff_fct==None:
            self.veff_fct=Function('spline',self.rgrid,self.veff)
        return self.veff_fct(r,der=der)
        
        
    def get_radial_density(self):
        return self.rgrid,self.dens        
        
        
    def Rnl(self,r,nl,der=0):
        """ Rnl(r,'2p') or Rnl(r,(2,1))"""
        nls=orbit_transform(nl,string=True)
        return self.Rnl_fct[nls](r,der=der)
        
        
    def unl(self,r,nl,der=0):
        """ unl(r,'2p')=Rnl(r,'2p')/r or unl(r,(2,1))..."""
        nls=orbit_transform(nl,string=True)
        return self.unl_fct[nls](r,der=der)
        
        
    def get_valence_orbitals(self):
        """ Get list of valence orbitals, e.g. ['2s','2p'] """
        return self.valence
    
    
    def get_symbol(self):
        """ Return atom's chemical symbol. """
        return self.symbol
        
        
    def get_comment(self):
        """ One-line comment, e.g. 'H, charge=0, frauenheim, r0=4' """
        comment='%s xc=%s charge=%.1f conf:%s' %(self.symbol,self.xc,float(self.charge),self.confinement_potential.get_comment())
        return comment
    
    
    def get_valence_energies(self):
        """ Return list of valence energies, e.g. ['2s','2p'] --> [-39.2134,-36.9412] """
        if not self.solved:
            raise AssertionError('run calculations first.')
        return [(nl,self.enl[nl]) for nl in self.valence]
    
    
    def write_functions(self,file,only_valence=True,step=1):
        """ Append functions (unl,v_effective,...) into file (only valence functions by default). 
        
        Parameters:
        -----------
        file: output file name (e.g. XX.elm)
        only_valence: output of only valence orbitals
        step: step size for output grid
        
        """
        if not self.solved:
            raise AssertionError('run calculations first.')
        if only_valence:
            orbitals=self.valence
        else:
            orbitals=[nl for n,l,nl in self.list_states()]
        o=open(file,'a')
        for nl in orbitals:
            print>>o, '\n\nu_%s=' %nl
            for r,u in zip(self.rgrid[::step],self.unlg[nl][::step]):
                print>>o, r,u
            
        print>>o,'\n\nv_effective='
        for r,ve in zip(self.rgrid[::step],self.veff[::step]):
                print>>o, r,ve        
        print>>o,'\n\nconfinement='
        for r,vc in zip(self.rgrid[::step],self.conf[::step]):
                print>>o, r,vc
        print>>o,'\n\n'
        
        
def shoot(u,dx,c2,c1,c0,N):       
    """
    Integrate diff equation
            
           2
         d u      du
         --- c  + -- c  + u c  = 0
           2  2   dx  1      0
         dx
         
    in equispaced grid (spacing dx) using simple finite difference formulas 
    
    u'(i)=(u(i+1)-u(i-1))/(2*dx) and
    u''(i)=(u(i+1)-2*u(i)+u(i-1))/dx**2
    
    u[0:2] *has already been set* according to boundary conditions.
    
    return u, number of nodes, the discontinuity of derivative at 
    classical turning point (ctp), and ctp
    c0(r) is negative with large r, and turns positive at ctp.
    """
    fp=c2/dx**2 + 0.5*c1/dx
    fm=c2/dx**2 - 0.5*c1/dx
    f0=c0-2*c2/dx**2
    
    # backward integration down to classical turning point ctp 
    # (or one point beyond to get derivative)                                
    # If no ctp, integrate half-way
    u[-1]=1.0
    u[-2]=u[-1]*f0[-1]/fm[-1]
    all_negative=all(c0<0)
    for i in xrange(N-2,0,-1): 
        u[i-1]=(-fp[i]*u[i+1]-f0[i]*u[i])/fm[i]
        if abs(u[i-1])>1E10: u[i-1:]*=1E-10 #numerical stability            
        if c0[i]>0: 
            ctp=i
            break
        if all_negative and i==N/2: 
            ctp=N/2
            break
                                   
    utp=u[ctp]
    utp1=u[ctp+1]
    dright=(u[ctp+1]-u[ctp-1])/(2*dx)
    
    for i in xrange(1,ctp+1):
        u[i+1]=(-f0[i]*u[i]-fm[i]*u[i-1])/fp[i]
    
    dleft=(u[ctp+1]-u[ctp-1])/(2*dx)
    scale=utp/u[ctp]    
    u[:ctp+1]*=scale
    u[ctp+1]=utp1 #above overrode
    dleft*=scale
    u=u*nu.sign(u[1])
        
    nodes=sum( (u[0:ctp-1]*u[1:ctp])<0 )
    A=(dright-dleft)*utp
    return u, nodes, A, ctp                                          
            
            
class RadialGrid:
    def __init__(self,grid):
        """ 
        mode
        ----
        
        rmin                                                        rmax
        r[0]     r[1]      r[2]            ...                     r[N-1] grid
        I----'----I----'----I----'----I----'----I----'----I----'----I
           r0[0]     r0[1]     r0[2]       ...              r0[N-2]       r0grid
           dV[0]     dV[1]     dV[2]       ...              dV[N-2]       dV
           
           dV[i] is volume element of shell between r[i] and r[i+1]
        """
        
        rmin, rmax=grid[0], grid[-1]
        N=len(grid)
        self.N=N
        self.grid=grid
        self.dr=self.grid[1:N]-self.grid[0:N-1]
        self.r0=self.grid[0:N-1]+self.dr/2
        # first dV is sphere (treat separately), others are shells 
        self.dV=4*nu.pi*self.r0**2*self.dr
        self.dV*=(4*nu.pi*rmax**3/3)/sum(self.dV)

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
        
    def get_comment(self):
        return self.comment        
        
    def none(self,r):
        return 0.0
        
    def frauenheim(self,r):
        return (r/self.r0)**2
        
    def __call__(self,r):
        return self.f(r)
        
    
        


class XC_PW92:
    def __init__(self):
        """ The Perdew-Wang 1992 LDA exchange-correlation functional. """
        self.small=1E-90
        self.a1 = 0.21370
        self.c0 = 0.031091
        self.c1 = 0.046644
        self.b1 = 1.0/2.0/self.c0*nu.exp(-self.c1/2.0/self.c0)
        self.b2 = 2*self.c0*self.b1**2
        self.b3 = 1.6382
        self.b4 = 0.49294
    
    def exc(self,n,der=0):
        """ Exchange-correlation with electron density n. """
        if n<self.small:
            return 0.0
        else:
            return self.e_x(n,der=der)+self.e_corr(n,der=der)
        
    def e_x(self,n,der=0):
        """ Exchange. """
        if der==0:
            return -3.0/4*(3*n/pi)**(1.0/3)
        elif der==1:
            return -3.0/(4*pi)*(3*n/pi)**(-2.0/3)

    def e_corr(self,n,der=0):
        """ Correlation energy. """
        rs = (3.0/(4*pi*n))**(1.0/3)
        aux=2*self.c0*( self.b1*sqrt(rs)+self.b2*rs+self.b3*rs**(3.0/2)+self.b4*rs**2 )
        if der==0:
            return -2*self.c0*(1+self.a1*rs)*log(1+aux**-1)
        elif der==1:
            return ( -2*self.c0*self.a1*log(1+aux**-1) \
                   -2*self.c0*(1+self.a1*rs)*(1+aux**-1)**-1*(-aux**-2)\
                   *2*self.c0*(self.b1/(2*sqrt(rs))+self.b2+3*self.b3*sqrt(rs)/2+2*self.b4*rs) )*( -(4*pi*n**2*rs**2)**-1 )
                   
    def vxc(self,n):
        """ Exchange-correlation potential (functional derivative of exc). """
        eps=1E-9*n       
        if n<self.small:
            return 0.0
        else: 
            return self.exc(n)+n*self.exc(n,der=1)

    
    
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
        

    

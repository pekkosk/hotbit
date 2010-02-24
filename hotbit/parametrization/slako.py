import numpy as nu
import pylab as pl
import sys
from box.timing import Timer
from util import tail_smoothening
from time import asctime
import math
sin=math.sin
cos=math.cos
tan=math.tan
arccos=math.acos
sqrt=math.sqrt
pi=math.pi
zeros=nu.zeros

class SlaterKosterTable:
    def __init__(self,ela,elb,txt=None,timing=False):
        """ Construct Slater-Koster table for given elements.
                
        parameters:
        -----------
        ela:    element objects (KSAllElectron or Element)
        elb:    element objects (KSAllElectron or Element)    
        txt:    output file object or file name
        timing: output of timing summary after calculation
        """
        self.ela=ela
        self.elb=elb
        self.timing=timing
        if txt==None:
            self.txt=sys.stdout
        else:
            if type(txt)==type(''):
                self.txt=open(txt,'a')
            else:
                self.txt=txt                
        self.comment=self.ela.get_comment()
        if ela.get_symbol()!=elb.get_symbol():
            self.nel=2
            self.pairs=[(ela,elb),(elb,ela)]
            self.elements=[ela,elb]
            self.comment+='\n'+self.elb.get_comment()
        else:
            self.nel=1
            self.pairs=[(ela,elb)]
            self.elements=[ela]
        self.timer=Timer('SlaterKosterTable',txt=self.txt,enabled=timing)
                                        
        print>>self.txt, '\n\n\n\n'                                        
        print>>self.txt, '************************************************'
        print>>self.txt, 'Slater-Koster table construction for %2s and %2s' %(ela.get_symbol(),elb.get_symbol())
        print>>self.txt, '************************************************'
        
    def __del__(self):
        self.timer.summary()          
        
        
    def get_table(self):
        """ Return tables. """
        return self.Rgrid, self.tables        
                
                
    def smooth_tails(self):
        """ Smooth the behaviour of tables near cutoff. """
        for p in range(self.nel):
            for i in range(20):
                self.tables[p][:,i]=tail_smoothening(self.Rgrid,self.tables[p][:,i])                
        
        
    def write(self,filename=None):
        """ Use symbol1_symbol2.par as default. """
        self.smooth_tails()
        if filename==None: fn='%s_%s.par' %(self.ela.get_symbol(),self.elb.get_symbol())
        else: fn=filename
        f=open(fn,'w')
        print>>f, 'slako_comment='
        print>>f, self.get_comment(), '\n\n'
        for p,(e1,e2) in enumerate(self.pairs):
            print>>f, '%s_%s_table=' %(e1.get_symbol(),e2.get_symbol())
            for i,R in enumerate(self.Rgrid):
                print>>f, '%.6e' %R,
                for t in xrange(20):
                    x=self.tables[p][i,t]
                    if abs(x)<1E-90:
                        print>>f, '0.',
                    else:
                        print>>f, '%.6e' %x,
                print>>f     
            print>>f, '\n\n'                                                                                   
        f.close()
    
    
    def plot(self,filename=None):
        """ Plot the Slater-Koster table with matplotlib. 
        
        parameters:
        ===========
        filename:     for graphics file
        
        """
        fig=pl.figure()
        fig.subplots_adjust(hspace=0.0001,wspace=0.0001)
        mx=max(1,self.tables[0].max())
        if self.nel==2:
            mx=max(mx,self.tables[1].max())
        for i in xrange(10):
            name=integrals[i]
            ax=pl.subplot(5,2,i+1)
            for p,(e1,e2) in enumerate(self.pairs):
                s1, s2=e1.get_symbol(), e2.get_symbol()
                if p==0: 
                    s='-'
                    lw = 1
                    alpha = 1.0
                else: 
                    s='--'
                    lw = 4
                    alpha = 0.2
                if all(abs(self.tables[p][:,i])<1E-10):
                    ax.text(0.03,0.02+p*0.15,'No %s integrals for <%s|%s>' %(name,s1,s2),transform=ax.transAxes,size=10)
                    pl.yticks([],[])
                    pl.xticks([],[])
                else:
                    pl.plot(self.Rgrid,self.tables[p][:,i],c='r',ls=s,lw=lw,alpha=alpha)
                    pl.plot(self.Rgrid,self.tables[p][:,i+10],c='b',ls=s,lw=lw,alpha=alpha)
                    pl.axhline(0,c='k',ls='--')
                    pl.title(name,position=(0.9,0.8)) 
                    if ax.is_last_row():
                        pl.xlabel('r (Bohr)')                                        
                    else:
                        pl.xticks([],[])
                    if not ax.is_first_col():                   
                        pl.yticks([],[])
                    pl.ylim(-mx,mx)
                    pl.xlim(0)
        
        pl.figtext(0.3,0.95,'H',color='r',size=20)
        pl.figtext(0.34,0.95,'S',color='b',size=20)
        pl.figtext(0.38,0.95,' Slater-Koster tables',size=20)
        e1, e2 =  self.ela.get_symbol(),self.elb.get_symbol()
        pl.figtext(0.3,0.92,'(thin solid: <%s|%s>, wide dashed: <%s|%s>)' %(e1,e2,e2,e1),size=10)
        
        file = '%s_%s_slako.pdf' %(e1,e2)
        if filename!=None:
            file = filename
        pl.savefig(file)            
    
    
    def get_comment(self):
        """ Get comments concerning parametrization. """        
        return self.comment
        
        
    def set_comment(self,comment):
        """ Add optional one-liner comment for documenting the parametrization. """
        self.comment+='\n'+comment
        
        
    def get_range(self,fractional_limit):
        """ Define ranges for the atoms: largest r such that Rnl(r)<limit. """
        self.timer.start('define ranges')
        wf_range=0.0
        for el in self.elements:
            r=max( [el.get_wf_range(nl,fractional_limit) for nl in el.get_valence_orbitals()] )
            print>>self.txt, 'wf range for %s=%10.5f' %(el.get_symbol(),r)
            wf_range=max(r,wf_range)
        if wf_range>20:
            raise AssertionError('Wave function range >20 Bohr. Decrease wflimit?')
        return wf_range
        self.timer.stop('define ranges')        
        
        
    def run(self,R1,R2,N,ntheta=150,nr=50,wflimit=1E-7):
        """ Calculate the Slater-Koster table. 
         
        parameters:
        ------------
        R1, R2, N: make table from R1 to R2 with N points
        ntheta: number of angular divisions in polar grid. (more dense towards bonding region)
        nr:     number of radial divisions in polar grid. (more dense towards origins)
                with p=q=2 (powers in polar grid) ntheta~3*nr is optimum (with fixed grid size)
                with ntheta=150, nr=50 you get~1E-4 accuracy for H-elements
                (beyond that, gain is slow with increasing grid size)
        wflimit: use max range for wfs such that at R(rmax)<wflimit*max(R(r))
        """
        if R1<1E-3:
            raise AssertionError('For stability; use R1>~1E-3')
        self.timer.start('calculate tables')   
        self.wf_range=self.get_range(wflimit)        
        Rgrid=nu.linspace(R1,R2,N)
        self.N=N
        self.Rgrid=Rgrid
        self.dH=0.0
        self.Hmax=0.0
        if self.nel==1: self.tables=[nu.zeros((N,20))]
        else: self.tables=[nu.zeros((N,20)),nu.zeros((N,20))]
        
        print>>self.txt, 'Start making table...'
        for Ri,R in enumerate(Rgrid):
            if R>2*self.wf_range: 
                break
            grid, areas = self.make_grid(R,nt=ntheta,nr=nr)
            if  Ri==N-1 or nu.mod(Ri,N/10)==0:                    
                    print>>self.txt, 'R=%8.2f, %i grid points ...' %(R,len(grid))
            for p,(e1,e2) in enumerate(self.pairs):
                selected=select_integrals(e1,e2) 
                if Ri==0:
                    print>>self.txt, 'R=%8.2f %s-%s, %i grid points, ' %(R,e1.get_symbol(),e2.get_symbol(),len(grid)),
                    print>>self.txt, 'integrals:', 
                    for s in selected: print>>self.txt, s[0],
                    print>>self.txt 
                
                S,H,H2=self.calculate_mels(selected,e1,e2,R,grid,areas)
                self.Hmax=max(self.Hmax,max(abs(H)))
                self.dH=max(self.dH,max(abs(H-H2)))
                self.tables[p][Ri,:10]=H
                self.tables[p][Ri,10:]=S
        
        print>>self.txt, 'Maximum value for H=%.2g' %self.Hmax
        print>>self.txt, 'Maximum error for H=%.2g' %self.dH        
        print>>self.txt, '     Relative error=%.2g %%' %(self.dH/self.Hmax*100)
        self.timer.stop('calculate tables')  
        self.comment+='\n'+asctime()
        self.txt.flush()
                    
                    
                    
                    
    def calculate_mels(self,selected,e1,e2,R,grid,area):
        """ 
        Perform integration for selected H and S integrals.
         
        parameters:
        -----------
        selected: list of [('dds','3d','4d'),(...)]
        e1: <bra| element
        e2: |ket> element
        R: e1 is at origin, e2 at z=R
        grid: list of grid points on (d,z)-plane
        area: d-z areas of the grid points.
        
        return:
        -------
        List of H,S and H2 for selected integrals. H2 is calculated using different
        technique and can be used for error estimation.
        
        S: simply R1*R2*angle-part
        H: operate (derivate) R2 <R1|t+Veff1+Veff2-Conf1-Conf2|R2>
        H2: operate with full h2 and hence use eigenvalue of |R2> with full Veff2
              <R1|(t1+Veff1)+Veff2-Conf1-Conf2|R2> 
            = <R1|h1+Veff2-Conf1-Conf2|R2> (operate with h1 on left)
            = <R1|e1+Veff2-Conf1-Conf2|R2> 
            = e1*S + <R1|Veff2-Conf1-Conf2|R2> 
            -> H and H2 can be compared and error estimated
        """
        self.timer.start('calculate_mels')
        Sl, Hl, H2l=nu.zeros(10), nu.zeros(10), nu.zeros(10)
        
        # common for all integrals (not wf-dependent parts)
        self.timer.start('prelude')
        N=len(grid)
        gphi, radii, v1, v2=zeros((N,10)), zeros((N,2)), zeros(N), zeros(N)
        for i,(d,z) in enumerate(grid):
            r1, r2=sqrt(d**2+z**2), sqrt(d**2+(R-z)**2)
            t1, t2=arccos(z/r1), arccos((z-R)/r2)
            radii[i,:]=[r1,r2]
            gphi[i,:]=g(t1,t2)
            v1[i]=e1.effective_potential(r1)-e1.confinement_potential(r1) 
            v2[i]=e2.effective_potential(r2)-e2.confinement_potential(r2) 
        self.timer.stop('prelude')                             
        
        # calculate all selected integrals
        for integral,nl1,nl2 in selected:           
            index=integrals.index(integral)
            S, H, H2=0.0, 0.0, 0.0
            l2=angular_momentum[nl2[1]]
            for i,dA in enumerate(area):            
                r1, r2=radii[i,:]
                d, z=grid[i]
                aux=gphi[i,index]*dA*d
                
                Rnl1, Rnl2, ddunl2=e1.Rnl(r1,nl1), e2.Rnl(r2,nl2), e2.unl(r2,nl2,der=2)
                
                S+= Rnl1*Rnl2*aux 
                H+= Rnl1*( -0.5*ddunl2/r2 + (v1[i]+v2[i]+l2*(l2+1)/(2*r2**2))*Rnl2 )*aux
                H2+= Rnl1*Rnl2*aux*( v2[i]-e1.confinement_potential(r1) )
            H2+=e1.get_epsilon(nl1)*S 
            Sl[index]=S
            Hl[index]=H
            H2l[index]=H2
            
        self.timer.stop('calculate_mels')
        return Sl, Hl, H2l
        
        
    def make_grid(self,Rz,nt,nr,p=2,q=2,view=False):
        """
        Construct a double-polar grid.
        
        Parameters:
        -----------
        Rz: element 1 is at origin, element 2 at z=Rz
        nt: number of theta grid points
        nr: number of radial grid points
        p: power describing the angular distribution of grid points (larger puts more weight 
           towards theta=0)
        q: power describing the radial disribution of grid points (larger puts more weight
           towards centers)   
        view: view the distribution of grid points with pylab.
          
        Plane at R/2 divides two polar grids.
                
                               
         ^ (z-axis)     
         |--------_____               phi_j
         |       /     ----__         *
         |      /            \       /  *              
         |     /               \    /  X *                X=coordinates of the center of area element(z,d), 
         |    /                  \  \-----* phi_(j+1)     area=(r_(i+1)**2-r_i**2)*(phi_(j+1)-phi_j)/2
         |   /                    \  r_i   r_(i+1)
         |  /                      \
         | /                       |
         *2------------------------|           polar centered on atom 2
         | \                       |
         |  \                     /                                                     1
         |   \                   /                                                     /  \
         |-------------------------- z=h -line         ordering of sector slice       /     \
         |   /                   \                                      points:      /        \
         |  /                     \                                                 /          \
         | /                       |                                               /     0       4
         *1------------------------|--->      polar centered on atom 1            2            /
         | \                       |    (r_perpendicular (xy-plane) = 'd-axis')    \        /
         |  \                      /                                                 \   /
         |   \                    /                                                    3
         |    \                  /
         |     \               /
         |      \           /
         |       \ ___ ---
         |---------
         
        """ 
        self.timer.start('make grid')
        rmin, rmax=(1E-7, self.wf_range)
        max_range=self.wf_range
        h=Rz/2
        T=nu.linspace(0,1,nt)**p*nu.pi
        R=rmin+nu.linspace(0,1,nr)**q*(rmax-rmin)
        
        grid=[]
        area=[]
        # first calculate grid for polar centered on atom 1:
        # the z=h-like starts cutting full elements starting from point (1)
        for j in xrange(nt-1):
            for i in xrange(nr-1):
                # corners of area element
                d1,z1=R[i+1]*sin(T[j]),R[i+1]*cos(T[j])
                d2,z2=R[i]*sin(T[j]),R[i]*cos(T[j])
                d3,z3=R[i]*sin(T[j+1]),R[i]*cos(T[j+1])
                d4,z4=R[i+1]*sin(T[j+1]),R[i+1]*cos(T[j+1])
                A0=(R[i+1]**2-R[i]**2)*(T[j+1]-T[j])/2
                
                if z1<=h:
                    # area fully inside region
                    r0=0.5*(R[i]+R[i+1])
                    t0=0.5*(T[j]+T[j+1])
                    A=A0
                elif z1>h and z2<=h and z4<=h:
                    # corner 1 outside region
                    Th=nu.arccos(h/R[i+1])
                    r0=0.5*(R[i]+R[i+1])
                    t0=0.5*(Th+T[j+1])
                    A=A0
                    A-=0.5*R[i+1]**2*(Th-T[j])-0.5*h**2*(tan(Th)-tan(T[j])) 
                elif z1>h and z2>h and z3<=h and z4<=h:
                    # corners 1 and 2 outside region
                    Th1=nu.arccos(h/R[i])
                    Th2=nu.arccos(h/R[i+1])
                    r0=0.5*(R[i]+R[i+1])
                    t0=0.5*(Th2+T[j+1])
                    A=A0
                    A-=A0*(Th1-T[j])/(T[j+1]-T[j])
                    A-=0.5*R[i+1]**2*(Th2-Th1)-0.5*h**2*(tan(Th2)-tan(Th1))
                elif z1>h and z2>h and z4>h and z3<=h:
                    # only corner 3 inside region
                    Th=nu.arccos(h/R[i])
                    r0=0.5*(R[i]+h/cos(T[j+1]))
                    t0=0.5*(Th+T[j+1])
                    A=0.5*h**2*(tan(T[j+1])-tan(Th)) - 0.5*R[i]**2*(T[j+1]-Th)
                elif z1>h and z4>h and z2<=h and z3<=h:
                    # corners 1 and 4 outside region
                    r0=0.5*(R[i]+h/cos(T[j+1]))
                    t0=0.5*(T[j]+T[j+1])
                    A=0.5*h**2*(tan(T[j+1])-tan(T[j])) - 0.5*R[i]**2*(T[j+1]-T[j])
                elif z3>h:
                    A=-1
                else:
                    raise RuntimeError('Illegal coordinates.')
                d,z=(r0*sin(t0),r0*cos(t0))
                if A>0 and sqrt(d**2+z**2)<max_range and sqrt(d**2+(Rz-z)**2)<max_range:
                    grid.append([d,z])
                    area.append(A)    
                                               
        self.timer.start('symmetrize')                                               
        # calculate the polar centered on atom 2 by mirroring the other grid                                               
        grid=nu.array(grid)
        area=nu.array(area)
        grid2=grid.copy()
        grid2[:,1]=-grid[:,1]
        shift=nu.zeros_like(grid)
        shift[:,1]=2*h
        grid=nu.concatenate((grid,grid2+shift))
        area=nu.concatenate((area,area))
        self.timer.stop('symmetrize')
                
        if view:
            pl.plot([h,h,h])        
            pl.scatter(grid[:,0],grid[:,1],s=10*area/max(area))
            pl.show()
            
        self.timer.stop('make grid')            
        return grid,area
        
        
        
        
integrals =['dds','ddp','ddd','pds','pdp','pps','ppp','sds','sps','sss']
angular_momentum={'s':0,'p':1,'d':2}


def select_orbitals(val1,val2,integral):
    """ 
    Select orbitals from given valences to calculate given integral. 
    e.g. ['2s','2p'],['4s','3d'],'sds' --> '2s' & '3d'
    """
    nl1=None
    for nl in val1:
        if nl[1]==integral[0]: nl1=nl
    nl2=None
    for nl in val2:
        if nl[1]==integral[1]: nl2=nl
    return nl1,nl2
        
        
def select_integrals(e1,e2):
    """ Return list of integrals (integral,nl1,nl2) to be done for element pair e1,e2. """
    selected=[]
    val1, val2 = e1.get_valence_orbitals(), e2.get_valence_orbitals()
    for ii,integral in enumerate(integrals):
        nl1, nl2=select_orbitals(val1,val2,integral)
        if nl1==None or nl2==None:
            continue
        else:
            selected.append( (integral,nl1,nl2) )
    return selected                
            
        
        
def g(t1,t2):
    """
    Return the angle-dependent part of the two-center 
    integral (it) with t1=theta_1 (atom at origin)
    and t2=theta2 (atom at z=Rz). These dependencies
    come after integrating analytically over phi.
    """
    c1, c2, s1, s2=cos(t1), cos(t2), sin(t1), sin(t2)   
    return nu.array([5.0/8*(3*c1**2-1)*(3*c2**2-1),\
            15.0/4*s1*c1*s2*c2,\
            15.0/16*s1**2*s2**2,\
            sqrt(15.0)/4*c1*(3*c2**2-1),\
            sqrt(45.0)/4*s1*s2*c2,\
            3.0/2*c1*c2,\
            3.0/4*s1*s2,\
            sqrt(5.0)/4*(3*c2**2-1),\
            sqrt(3.0)/2*c2,\
            0.5])
           
        
if __name__=='__main__':
    print select_orbital2(['2s','2p'],['4s','3d'],'sss')
    

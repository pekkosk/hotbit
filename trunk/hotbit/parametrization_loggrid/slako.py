import numpy as nu
import pylab as pl
import sys
sin=nu.sin
cos=nu.cos
tan=nu.tan
sqrt=nu.sqrt
pi=nu.pi

integrals =['dds','ddp','ddd','pds','pdp','pps','ppp','sds','sps','sss']
integral_l=[(2,2),(2,2),(2,2),(1,2),(1,2),(1,1),(1,1),(0,2),(0,1),(0,0)]
angular_momentum={'s':0,'p':1,'d':2}

def select_orbitals(val1,val2,integral):
    """ 
    Select orbitals from given valences to calculate given integral. 
    e.g. ['2s','2p'],['4s','3d'],'sds' --> '2s' & '3d'
    """
    nl1=None
    for nl in val1:
        if nl[1]==integral[0]:
            nl1=nl
    nl2=None
    for nl in val2:
        if nl[1]==integral[1]:
            nl2=nl
    return nl1,nl2
        
        
        
def g(i,t1,t2):
    """
    Return the angle-dependent part of the two-center 
    integral (it) with t1=theta_1 (atom at origin)
    and t2=theta2 (atom at z=Rz). These dependencies
    come after integrating analytically over phi.
    """
    if i=='dds':
        return 5.0/8*(3*cos(t1)**2-1)*(3*cos(t2)**2-1)
    elif i=='ddp':
        return 15.0/4*sin(t1)*cos(t1)*sin(t2)*cos(t2)
    elif i=='ddd':
        return 15.0/16*sin(t1)**2*sin(t2)**2
    elif i=='pds':
        return sqrt(15.0)/4*cos(t1)*(3*cos(t2)**2-1)
    elif i=='pdp':
        return sqrt(45.0)/4*sin(t1)*sin(t2)*cos(t2)
    elif i=='pps':
        return 3.0/2*cos(t1)*cos(t2)
    elif i=='ppp':
        return 3.0/4*sin(t1)*sin(t2)
    elif i=='sds':
        return sqrt(5.0)/4*(3*cos(t2)**2-1)
    elif i=='sps':
        return sqrt(3.0)/2*cos(t2)
    elif i=='sss': 
        return 0.5
 
    

class SlaterKosterTable:
    def __init__(self,ela,elb,limit=1E-10,out=None):
        """ Given atoms as Element objects. """
        self.ela=ela
        self.elb=elb
        if out==None:
            self.out=sys.stdout
        else:
            self.out=out
        self.comment=self.ela.get_comments()+'\n'+self.elb.get_comments()
        if ela.get_symbol()!=elb.get_symbol():
            self.nel=2
            self.pairs=[(self.ela,self.elb),(self.elb,self.ela)]
            self.elements=[self.ela,self.elb]
        else:
            self.nel=1
            self.pairs=[(self.ela,self.elb)]
            self.elements=[self.ela]
        self.define_ranges(limit)
        
    def write_par(self,name=None):
        """ Use symbol1_symbol2.par as default. """
        pass
    
    def plot(self):
        """ Plot table. """
        mx=self.tables[0].max()
        if self.nel==2:
            mx=max(mx,self.tables[1].max())
        for i in range(10):
            pl.subplot(5,2,i+1)
            pl.plot(self.Rgrid,self.tables[0][:,i],label='H')
            pl.plot(self.Rgrid,self.tables[0][:,i+10],label='S')
            if self.nel==2:
                pl.plot(self.Rgrid,self.tables[1][:,i],label='H')
                pl.plot(self.Rgrid,self.tables[1][:,i+10],label='S')
            pl.ylim(-mx,mx)
            pl.xlim(0)
            pl.title('table')
            pl.xlabel('r (Bohr)')
            pl.ylabel('H (Ha) and S')
        pl.show()
    
    def get_comment(self):
        """ Get comments concerning parametrization. """
        return self.comment
        
    def set_comment(self,comment):
        """ Add optional one-liner comment for documenting the parametrization. """
        self.comment+='\n'+comment
        
    def define_ranges(self,limit):
        """ Define ranges for the atoms: largest r such that Rnl(r)<limit. """
        grid=nu.linspace(30,0,300)
        self.rnl_range=0.0
        for el in self.elements:
            for r in grid:
                if max( [abs(el.Rnl(r,nl)) for nl in el.get_valence_orbitals()] )>limit:
                    el.Rnl_range=r
                    self.rnl_range=max(r,self.rnl_range)
                    print>>self.out, 'Rnl range for %s=%10.5f' %(el.get_symbol(),r)
                    break
        
        
    def calculate_tables(self,R1,R2,N,ntheta=100,nr=50):
        """ 
        Calculate the whole table from R1 to R2 with N points.
        rmin is closest r to nuclei and rmax is the maximum radius in the double-polar grid.
        ntheta and nr are the number of angular and radial grid points in double-polar grid.
        """
        
        Rgrid=nu.linspace(R1,R2,N)
        self.Rgrid=Rgrid
        self.dH=0.0
        self.Hmax=0.0
        self.tables=[nu.zeros((N,20))]*self.nel
        for Ri,R in enumerate(Rgrid):
            grid,areas=self.make_grid(R,rmin=1E-7,rmax=self.rnl_range,nt=ntheta,nr=nr,max_range=self.rnl_range)
            for p,(e1,e2) in enumerate(self.pairs):
                val1=e1.get_valence_orbitals()
                val2=e2.get_valence_orbitals()
                for ii,integral in enumerate(integrals):
                    # choose which orbitals to use for this integral
                    nl1,nl2=select_orbitals(val1,val2,integral)
                    if nl1==None or nl2==None:
                        continue
                    print 'R=%8.2f <%s(%s)|..|%s(%s)> %s-integral, %i grid points' \
                           %(R,e1.get_symbol(),nl1,e2.get_symbol(),nl2,integral,len(grid))
                    S,H,H2=self.integrate(integral,e1,nl1,e2,nl2,R,grid,areas)
                    self.Hmax=max(self.Hmax,abs(H))
                    self.dH=max(self.dH,abs(H-H2))
                    self.tables[p][Ri,ii],self.tables[p][Ri,ii+10]=H,S
        print>>self.out, 'Maximum uncertainty for H-element=%.2g' %self.dH
        print>>self.out, 'Maximum value of |<H>|=%.2g' %self.Hmax
        print>>self.out, 'Relative uncertainty=%.2g' %(self.dH/self.Hmax)
                    
                    
    def integrate(self,integral,e1,nl1,e2,nl2,R,grid,area):
        """ 
        Perform integration with nl1 orbital of element1 and nl2 orbital of element2,
        with element1 at origin and element2 ar z=R. Both H ans S elements.
        """
        S,H,H2=0.0,0.0,0.0
        l2=angular_momentum[nl2[1]]
        for (d,z),dA in zip(grid,area):
            r1=sqrt(d**2+z**2)
            r2=sqrt(d**2+(R-z)**2)
            theta1=nu.arccos(z/r1)
            theta2=nu.arccos((z-R)/r2)
            aux=g(integral,theta1,theta2)*dA*d
            v1, v2=e1.v_effective(r1)-e1.confinement(r1), e2.v_effective(r2)-e2.confinement(r2)
            Rnl1, Rnl2, ddunl2=e1.Rnl(r1,nl1), e2.Rnl(r2,nl2), e2.unl(r2,nl2,der=2)
            S+= Rnl1*Rnl2*aux
            H+= Rnl1*( -0.5*ddunl2/r2+(v1+v2+l2*(l2+1)/(2*r2**2))*Rnl2 )*aux
            H2+= Rnl1*Rnl2*aux*( e2.v_effective(r2)-e1.confinement(r1)-e1.confinement(r2) )
        H2+=e1.get_epsilon(nl1)*S 
        return S,H,H2
        
        
    def make_grid(self,Rz,rmin,rmax,nt=10,nr=10,p=2,q=2,max_range=100,view=False):
        """
        Construct a double-polar grid.
        
        Center 1 is at origin, center 2 at z=R. Plane at R/2 divides two polar grids.
        * nt=number of theta grid points, nr=number of radial grid points
        * p=power describing the angular distribution of grid points (larger puts more weight 
          towards theta=0)
        * q=power describing the radial disribution of grid points (larger puts more weight
          towards center)
        * max_range is the maximum range of grid point from either center.
        
                               
         ^ (z-axis)     
         |--------_____               phi_j
         |       /     ----__         *
         |      /            \       /  *              
         |     /               \    /  X *                X=coordinates of the center of area element(z,d), 
         |    /                  \  \-----* phi_(j+1)     area=(r_(i+1)**2-r_i**2)*(phi_(j+1)-phi_j)/2
         |   /                    \  r_i   r_(i+1)
         |  /                      \
         | /                       |
         *1------------------------|           polar centered on atom 2
         | \                       |
         |  \                     /                                                     1
         |   \                   /                                                     /  \
         |-------------------------- z=h -line                                        /     \
         |   /                   \                                                   /        \
         |  /                     \                                                 /          \
         | /                       |                                               /     0       4
         *1------------------------|--->      polar centered on atom 1            2            /
         | \                       |    (r_perpendicular (xy-plane) = d-axis)      \        /
         |  \                      /                                                 \   /
         |   \                    /                                                    3
         |    \                  /
         |     \               /
         |      \           /
         |       \ ___ ---
         |---------
         
        """ 
        h=Rz/2
        T=nu.linspace(0,1,nt)**p*nu.pi
        R=rmin+nu.linspace(0,1,nr)**q*(rmax-rmin)
        
        grid=[]
        area=[]
        for j in range(nt-1):
            for i in range(nr-1):
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
                                               
        grid=nu.array(grid)
        area=nu.array(area)
        grid2=grid.copy()
        grid2[:,1]=-grid[:,1]
        shift=nu.zeros_like(grid)
        shift[:,1]=2*h
        grid=nu.concatenate((grid,grid2+shift))
        area=nu.concatenate((area,area))
        
        #assert abs(pi*rmax**2-rmax**2*nu.arccos(h/rmax)+h*nu.sqrt(rmax**2-h**2)-sum(area))<1E-9
                
        if view:
            pl.plot([h,h,h])        
            pl.scatter(grid[:,0],grid[:,1],s=10*area/max(area))
            pl.show()
        return grid,area
        
        
if __name__=='__main__':
    print select_orbital2(['2s','2p'],['4s','3d'],'sss')
    
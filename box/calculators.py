import numpy as npy
import pylab as pl
pi=npy.pi
vec=npy.array

def V_LEPS(rab,rbc):
    """ 
    One model potential for NEB calculations that looks like this:
    
    .A...||||||........
    ......||||.........
    .......||..........
    ...................
    .......|...........
    .......||..........
    ......||||.........
    .....||||||......B.
    ....||||||||.......
    """
    def Q(r,alpha,r0,d): return d/2.*( 3./2*npy.exp(-2*alpha*(r-r0))-npy.exp(-alpha*(r-r0)) )
    def J(r,alpha,r0,d): return d/4.*( npy.exp(-2*alpha*(r-r0))-6*npy.exp(-alpha*(r-r0)) )
    
    a=0.05
    b=0.3
    c=0.05
    dab=4.746
    dbc=4.746
    dac=3.445
    r0=0.742
    alpha=1.942
    ai=1.0/(1+a)
    bi=1.0/(1+b)
    ci=1.0/(1+c)
    Qab=Q(rab,alpha,r0,dab)
    Qbc=Q(rbc,alpha,r0,dbc)
    Qac=Q(rab+rbc,alpha,r0,dac)
    Jab=J(rab,alpha,r0,dab)
    Jbc=J(rbc,alpha,r0,dbc)
    Jac=J(rab+rbc,alpha,r0,dac)
    return Qab*ai+Qbc*bi+Qac*ci-npy.sqrt(Jab**2*ai**2+Jbc**2*bi**2+Jac**2*ci**2\
                                        -Jab*Jbc*ai*bi-Jbc*Jac*bi*ci-Jab*Jac*ai*ci)
                                                
def neb_model_function(rab,x):                           
    rac=3.742
    kc=0.2025
    c=1.154
    return V_LEPS(rab,rac-rab)+2*kc*(rab-(rac/2-x/c))**2


def neb_model_1(x,y,der=0):
    f=neb_model_function
    d=1E-6
    if der==0:
        return f(x,y)
    elif der==1:
        return [(f(x+d,y)-f(x,y))/d,(f(x,y+d)-f(x,y))/d]        
        
    
    
    

class Calculator:
    def __init__(self,potential):
        
        if potential=='neb_model_1':
            self.f=neb_model_1
        self.it=0
        
    def __call__(self,r):
        return self.f(r)
    
    def get_potential_energy(self,atoms):
        r=atoms.get_positions()
        return self.f(r[0,0],r[0,1])
        
    def get_forces(self,atoms):
        r=atoms.get_positions()
        der=self.f(r[0,0],r[0,1],der=1)
        f0=-vec([der[0],der[1],0.0])
        if len(r)==1:
            return vec([f0])
        else:
            return vec([f0,vec([0]*3)*(len(r)-1)])
        
    def get_stress(self,atoms):
        return None
    
    def plot(self,extent,scatter=None,plot=None,out='screen'):
        x=npy.linspace(extent[0],extent[1],100)
        y=npy.linspace(extent[2],extent[3],100)
        X,Y=npy.meshgrid(x,y)
        Z=self.f(X,Y)
        im = pl.imshow(Z, interpolation='bilinear', origin='lower',
                       cmap=pl.cm.gray, extent=extent)
    
        minz,maxz=(vec(Z).min(),vec(Z).max())
        cset = pl.contour(Z, npy.arange(minz,maxz,(maxz-minz)/20),
                          origin='lower',linewidths=1,extent=extent)
        
        if scatter!=None:
            pl.scatter(scatter[0,:],scatter[1,:],color='r')
        if plot!=None:
            pl.plot(plot[0,:],plot[1,:],color='y')
            pl.scatter(plot[0,:],plot[1,:],color='y')
            
        if out=='screen':
            pl.show()
        else:
            pl.savefig(out)
            pl.clf()
        self.it+=1

        
if __name__=='__main__':
    calc=Calculator('neb_model_1')
    calc.plot((0.5,3.2,-2.0,1.5))
      
    
    
    
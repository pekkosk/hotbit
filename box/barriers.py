import os
from box import Atoms
from box import md
import numpy as nu
import ase
from box.interpolation import VectorSplineFunction
#from box.calculators import Calculator
from scipy.linalg import norm
import pylab as pl
vec=nu.array


class LEPS_Calculator:
    """ ase Calculator class """
    
    def __init__(self):
        """
        Construct calculator for ase.Atoms class using
        V_LEPS potential from NEB papers. It looks like this:
        
        .A...||||||........
        ......||||.........
        .......||..........
        ...................
        .......|...........
        .......||..........
        ......||||.........
        .....||||||......B.
        ....||||||||.......
        
        This is designed only for a single atom -Atoms-object.
        """        
        self.it=0
           
    def _V_LEPS(self,rab,rbc):
        """ 
        The complicated part in the function definition.
        """
        def Q(r,alpha,r0,d): return d/2.*( 3./2*nu.exp(-2*alpha*(r-r0))-nu.exp(-alpha*(r-r0)) )
        def J(r,alpha,r0,d): return d/4.*( nu.exp(-2*alpha*(r-r0))-6*nu.exp(-alpha*(r-r0)) )
        
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
        return Qab*ai+Qbc*bi+Qac*ci-nu.sqrt(Jab**2*ai**2+Jbc**2*bi**2+Jac**2*ci**2\
                                            -Jab*Jbc*ai*bi-Jbc*Jac*bi*ci-Jab*Jac*ai*ci)
                                                    
    def aux(self,rab,xx): 
        rac=3.742
        kc=0.2025
        c=1.154
        return self._V_LEPS(rab,rac-rab)+2*kc*(rab-(rac/2-xx/c))**2
                                                            
    def _function(self,x,y,der=0):
        """ LEPS potential written for 'atomic' positions x and y. """        
        # use auxiliary function
        d=1E-6
        if der==0:
            return self.aux(x,y)
        elif der==1:
            return [(self.aux(x+d,y)-self.aux(x,y))/d,(self.aux(x,y+d)-self.aux(x,y))/d]                    
 
    def get_potential_energy(self,atoms):
        r=atoms.get_positions()
        return self._function(r[0,0],r[0,1])
 
    def __call__(self,r):
        return self.f(r)

    def set_atoms(self,atoms):
        assert len(atoms)==1       
        
    def get_forces(self,atoms):
        """ Forces from LEPS potential for the one atom. """
        r=atoms.get_positions()
        der=self._function(r[0,0],r[0,1],der=1)
        f0=-vec([der[0],der[1],0.0])
        return vec([f0])
        #if len(r)==1:
            #return vec([f0])
        #else:
            #return vec([f0,vec([0]*3)*(len(r)-1)])
        
    def get_stress(self,atoms):
        """ This just has to exist. """
        return None
    
    def plot(self,scatter=None,screen=True,out=None):
        """ Plot the VLEP landscape. """
        X = nu.linspace(0.5,3.2,100)
        Y = nu.linspace(-2,1.4,100)
        Z=nu.zeros((100,100))
        for i,x in enumerate(X):
            for j,y in enumerate(Y):
                Z[j,i]=self._function(x,y)
                
        pl.contourf(X,Y,Z,100)
        pl.hot()
        
        if scatter!=None:
            pl.scatter(scatter[0,:],scatter[1,:],color='b')
        #if plot!=None:
            #pl.plot(plot[0,:],plot[1,:],color='y')
            #pl.scatter(plot[0,:],plot[1,:],color='y')
            
        if screen==True:
            pl.show()
        else:
            assert out!=None
            pl.savefig(out)
            pl.clf()
        self.it+=1
        



class MEP:
    def __init__(self,images,calc):
        self.images=images
        self.N=len(images[0])
        self.M=len(images)  
        self.dim=self.M*self.N*3     
        self.calc=calc 
        self.it=0
        self.dir='MEP'
        if not os.path.isdir(self.dir):
            os.system('mkdir %s' %self.dir)
        self.lock=False   
        self.lambdas=None
                 
                 
    def has(self,key):
        return False                 
        
        
    def __len__(self):
        return self.M*self.N
    
    #def set_optimizer(self,optimizer):
        #""" Set the optimizer used to optimize the whole path. """
        #self.optimizer=optimizer
    
    def interpolate(self):
        """ Use linear interpolation for the images. """
        first=self.images[0].get_positions()
        last=self.images[-1].get_positions()
        for i in range(1,self.M-1):
            self.images[i].set_positions(first+i*1.0/(self.M-1.0)*(last-first))                
                                          
                     
    def lock_positions(self):
        self.lock=True
                             
                             
    def release_positions(self):
        self.lock=False
                                     
    
    def set_positions(self,positions):
        """ Set positions for images, given positions for M*N 'atoms'. """
        if self.lock: 
            self.positions2=nu.zeros((self.M,self.N,3))
            n1 = 0
            for i in range(len(self.images)):
                n2 = n1 + self.N
                self.positions2[i] = positions[n1:n2]
                n1 = n2
            self.positions2.shape=(self.M,self.N*3)
            return
            
        n1 = 0
        for image in self.images:
            n2 = n1 + self.N
            image.set_positions(positions[n1:n2])
            n1 = n2
            
    def get_lambda_guess(self):
        #if self.lambdas==None:
        self.it+=1
        R=self.get_positions()
        #R.shape=(self.M,self.N,3)
        out='%s/landscape_%04i.png' %(self.dir,self.it)
        #print R
        #R.shape=(self.M,self.N,3)
        #print R
        #assert False
        #print vec([R[:,0],R[:,1]])
        
        self.calc.plot(scatter=vec([R[:,0],R[:,1]]),screen=False,out=out)
        
        if self.lambdas==None:
            return [0.]*(self.M-1)
        else:
            return self.lambdas      
                  
    def get_max_constraint_deviation(self):
        """ Return max( |R[i]-R[i-1]|-L/(M-1) ) """
        R=self.get_positions()
        R.shape=(self.M,self.N*3)     
        
        D=nu.zeros(self.M-1)
        for i in range(self.M-1):
            D[i] = nu.linalg.norm(R[i+1]-R[i]) 
        L1=sum(D)/(self.M-1)            
        return max( abs(D-L1) )/L1
        
                          
            
    def get_constraint_deviations(self):
        """ Return how much we deviatie from the M-1 constraints. """
        dev=[]
        L=0.0
        try:
            R=self.positions2
        except:
            R=self.get_positions()
            R.shape=(self.M,self.N*3)      
        for i in range(self.M-1):
            Di=nu.linalg.norm(R[i+1]-R[i])
            L+=Di
            dev.append( Di )
        return nu.array(dev) - L/(self.M-1)           
        
        
    def get_positions(self):
        """ Get path positions, as for M*N 'atoms'. """
        positions = nu.empty((self.M*self.N,3))
        n1 = 0
        for image in self.images:
            n2 = n1 + self.N
            positions[n1:n2] = image.get_positions()
            n1 = n2
        return positions           
        
        
    def get_potential_energy(self):
        """ Return value of the objective function (sum of all energies). """
        return sum( [self.calc.get_potential_energy(image) for image in self.images] )
        
        
    def _get_potential_forces(self):
        """ Return forces only from potential energy. 
        (first and last images should be optimized already.)
        """
        F=[vec([0.0]*3*self.N)]         
        for image in self.images[1:-1]:
            F.append(self.calc.get_forces(image).reshape(3*self.N))            
        F.append(vec([0.0]*3*self.N)) 
        return vec(F)
        
        
    def _get_constraint_forces(self,lambdas):
        """ Return the constraint forces for given lambdas (Lagrange multipliers). """ 
        assert len(lambdas)==self.M-1               
        R=self.get_positions()
        R.shape=(self.M,self.N*3)     
        
        #M=nu.zeros((self.M-1,self.M-1))
        #for i in range(self.M-1):
            #for j in range(self.M-1):
                #M[i,j]=
        
        d=[]
        for i in range(self.M-1):
            D=R[i+1]-R[i]
            d.append( D/nu.linalg.norm(D) )           
               
        lambda0=sum(lambdas)/(self.M-1)               
        F=[vec([0.0]*3*self.N)]         
        for k in range(1,self.M-1):
            F.append( d[k-1]*(lambdas[k-1]-lambda0) - d[k]*(lambdas[k]-lambda0) )
        F.append(vec([0.0]*3*self.N)) 
        return vec(F)
                    
    def set_lambdas(self,lambdas):
        self.lambdas=lambdas                    
                    
    def get_forces(self):
        f0=self._get_potential_forces()
        fc=self._get_constraint_forces(self.lambdas)
        F = f0 + fc
        F.shape=(self.M*self.N,3)   
        return F
                            
        
        
        
        
                
        #k=1000
        #for i in range(1,self.M-1):
            #tangent=ip.normalized_tangent(t[i])
            #F_homo=-k*nu.vdot((R[i]-Rh[i]),tangent)*tangent
            #F_tan=nu.vdot(tangent,F[i])*tangent
            #print 'homo:',norm(F_homo),'tan:',norm(F_tan)
            #F[i]=F[i]-F_tan+F_homo
            
        #F=vec(F)
        #F.shape=(self.M*self.N,3)
        #self.it+=1
        #rr=ip.homogenize(200)
        #Rh.shape=(self.M,self.N*3)
        #out='%s/landscape_%04i.png' %(self.dir,self.it)
        #self.calc.plot((0.5,3.2,-2.0,1.5),scatter=vec([R[:,0],R[:,1]]),\
                        #plot=vec([Rh[:,0],Rh[:,1]]),out=out)
        #ip.plot(out='%s/interpol_%04i.png' %(self.dir,self.it))
        
        #return vec(F)
             

        
        
        

class BTI:
    def __init__(self,images,calc):
        self.images=images
        self.N=len(images[0])
        self.M=len(images)  
        self.dim=self.M*self.N*3     
        self.calc=calc 
        self.it=0
        self.dir='MEP'
        if not os.path.isdir(self.dir):
            os.system('mkdir %s' %self.dir)
        
    def __len__(self):
        return self.M*self.N
        
    def get_energy(self,image):
        """ Return energy for given image. """
        
    def initial_interpolation(self):
        first=self.images[0].get_positions()
        last=self.images[-1].get_positions()
        for i in range(1,self.M-1):
            self.images[i].set_positions(first+i*1.0/(self.M-1.0)*(last-first))    
                 
          
    def get_forces(self):
        R=self.get_positions()
        R.shape=(self.M,self.N*3)
        ip=VectorSplineFunction(R,k=3,s=0)
        Rh=ip.homogenize()
        t=ip.closest_parameters(Rh)
        
        # draw tangents
        #if False:
        tangents=vec(ip.tangents())
        tangents.shape=(self.M,self.N,3)
        for image,i in zip(self.images,range(self.M)):
            image.set('tangent',tangents[i,:,:])
            image.write_vtk('%s/image_%04i_%04i.vtk' %(self.dir,self.it,i))
            
        F=[vec([0.0]*3*self.N)]
        for image in self.images[1:-1]:
            image.set_calculator(self.calc)
            F.append(image.get_forces().reshape(3*self.N))            
        F.append(vec([0.0]*3*self.N))
        
        k=1000
        for i in range(1,self.M-1):
            tangent=ip.normalized_tangent(t[i])
            F_homo=-k*nu.vdot((R[i]-Rh[i]),tangent)*tangent
            F_tan=nu.vdot(tangent,F[i])*tangent
            print 'homo:',norm(F_homo),'tan:',norm(F_tan)
            F[i]=F[i]-F_tan+F_homo
            
        F=vec(F)
        F.shape=(self.M*self.N,3)
        self.it+=1
        rr=ip.homogenize(200)
        Rh.shape=(self.M,self.N*3)
        out='%s/landscape_%04i.png' %(self.dir,self.it)
        self.calc.plot((0.5,3.2,-2.0,1.5),scatter=vec([R[:,0],R[:,1]]),\
                        plot=vec([Rh[:,0],Rh[:,1]]),out=out)
        ip.plot(out='%s/interpol_%04i.png' %(self.dir,self.it))
        
        return vec(F)
             
                
    def get_potential_energy(self):
        E=[]
        for image in self.images:
            image.set_calculator(self.calc)
            E.append(image.get_potential_energy())
        return max(E)
       
    def set_positions(self,positions):
        n1 = 0
        for image in self.images:
            n2 = n1 + self.N
            image.set_positions(positions[n1:n2])
            n1 = n2
        
    def get_positions(self):
        positions = nu.empty((self.M*self.N,3))
        n1 = 0
        for image in self.images:
            n2 = n1 + self.N
            positions[n1:n2] = image.get_positions()
            n1 = n2
        return positions   
        
        
        
if __name__=='__main__':
    fmax=1E-4
    calc=Calculator(potential='neb_model_1')
    first=Atoms(symbols='H',positions=[(3.0,-1.3,0.0)]) #,(0,0,0)])
    first.set_calculator(calc)
    md.quench_atoms(first,fmax=fmax)
    print 'first minimum:',first.get_positions()[0]
    
    
    
    last=Atoms('H',[(0.77,1.31,0.0)]) #,(0,0,0)])
    last.set_calculator(calc)
    md.quench_atoms(last,fmax=fmax)
    print 'last minimum:',last.get_positions()[0]
    
    
    images=[first.copy() for x in range(10)]
    images.append(last)
    
    
    bti=BTI(images,calc)
    bti.initial_interpolation()
    qn=ase.MDMin(bti,dt=0.05)
    writer=md.TrajectoryWriter(images,name='test')
    qn.attach(writer)
    qn.run(fmax=fmax)

    
    
    
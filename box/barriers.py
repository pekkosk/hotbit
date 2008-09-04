import os
from box import Atoms
from box import md
import numpy as npy
import ase
from box.interpolation import VectorSplineFunction
from box.calculators import Calculator
from scipy.linalg import norm
vec=npy.array

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
            F_homo=-k*npy.vdot((R[i]-Rh[i]),tangent)*tangent
            F_tan=npy.vdot(tangent,F[i])*tangent
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
        positions = npy.empty((self.M*self.N,3))
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

    
    
    
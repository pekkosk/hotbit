from ase import Atoms as ase_Atoms
import numpy as nu
from box.mix import  phival
from math import sin,cos


class WedgeAtoms(ase_Atoms):     

    def __init__(self, *args, **kwargs):
        ase_Atoms.__init__(self, *args, **kwargs)
    
               
    def set(self,omega):
        # TODO: disable set_pbc and set_cell
        self.omega=omega
        self.pbc=nu.array([True,False,False],bool)
        M1 = round(2*nu.pi/omega)
        if nu.abs(M1-2*nu.pi/omega)>1E-12 or omega>nu.pi:
            raise AssertionError('omega not 2*pi/M !')
        self.M1=int(M1)
        i = M1/2
        if nu.mod(M1,2)==1:
            self.ranges = nu.array([[-i,i],[0,0],[0,0]],int)
        else:
            self.ranges = nu.array([[-i+1,i],[0,0],[0,0]],int)
        
    
    def get_ranges(self):
        '''
        Return the ranges for symmetry operations.
        '''
        return self.ranges.copy()        
    
    
    def _check(self,n):
        '''
        Check that the symmetry operation is allowed.
        
        @param n: 3-tuple for transformation
        '''
        for i in range(3):
            if not self.ranges[i,0]<=n[i]<=self.ranges[i,1]:
                raise AssertionError('Illegal symmetry operation: %i %i %i' %n)
        
    
    def transform(self,r,n):
        '''
        Transform position r with S(n)
        
        @param r: position vector
        @param n: 3-tuple for symmetry operations
        '''
        self._check(n)
        rn=nu.zeros_like(r)
        x, y = r[0], r[1]
        rad = nu.sqrt(x**2+y**2)
        phi = phival(x,y)
        rn[0] = rad * cos( phi+n[0]*self.omega )
        rn[1] = rad * sin( phi+n[0]*self.omega )
        rn[2] = r[2]
        return rn
        
    
    def dtensor(self,r,n):
        '''
        Return the dyadic tensor 
        
                    d (R_j^n)_a
        T(jn)_ab = ------------- hat_a hat_b 
                    d (R_j)b
                    
        @param r: 
        @param n:
        '''
        self._check(n)
        rn = self.transform(r,n)
        rad = nu.sqrt(r[0]**2+r[1]**2)
        assert rad>1E-10
        T = nu.array([[r[0]*rn[0]+r[1]*rn[1],-(r[0]*rn[1]-r[1]*rn[0]),0],\
                      [r[0]*rn[1]-r[1]*rn[0],  r[0]*rn[0]+r[1]*rn[1] ,0],\
                      [       0,                      0,         rad**2]])/rad**2
        return T
        
    
    def axis_rotation(self,n):
        '''
        Return the rotation matrix for given symmetry operation.
        
        @param n: 3-tuple
        '''
        self._check(n)
        angle = n[0]*self.omega
        R = nu.array([[cos(angle),sin(angle),0],[-sin(angle),cos(angle),0],[0,0,1]])
        return R     
        
        
    def __eq__(self,other):
        if isinstance(other,WedgeAtoms):
            # More strict comparison for BravaisAtoms
            # TODO: check additional stuff for other generalized classes
            # (angles, torsions ...)
            if (self.positions-other.positions<1E-13).all() and \
               (self.pbc==other.pbc).all() and \
               self.get_chemical_symbols()==other.get_chemical_symbols() and \
               (self.cell==other.cell).all():
                return True
            else:
                return False
        else:
            # other is probably normal ase.Atoms; more loose check
            if (self.positions-other.positions<1E-13).all() and \
               (self.pbc==other.pbc).all() and \
               self.get_chemical_symbols()==other.get_chemical_symbols() and \
               (self.cell==other.cell).all():
                return True
            else:
                return False






class BravaisAtoms(ase_Atoms):     

    def __init__(self, *args, **kwargs):
        ase_Atoms.__init__(self, *args, **kwargs)
        ranges = []
        for i in range(3):
            if self.pbc[i]:
                ranges.append([-nu.Inf,nu.Inf])
            else:
                ranges.append([0,0])
        self.ranges = nu.array(ranges)
               
    
    def transform(self,r,n):
        '''
        Transform position r with S(n)
        
        @param r: position vector
        @param n: 3-tuple for symmetry operations
        '''
        # TODO: check than n belongs to allowed transformations
        rn=r.copy()
        for a in range(3):
            rn = rn + n[a]*self.cell[a,:]
        return rn
    
    
    def get_ranges(self):
        '''
        Return the ranges for symmetry operations.
        '''
        return self.ranges.copy() 
    
    
    def dtensor(self,r,n):
        '''
        Return the dyadic tensor 
        
                    d (R_j^n)_a
        T(jn)_ab = ------------- hat_a hat_b 
                    d (R_j)b
                    
        @param r: 
        @param n:
        '''
        return nu.eye(3)
    
    
    def axis_rotation(self,n):
        '''
        Return the rotation matrix for given symmetry operation.
        
        @param n: 3-tuple
        '''
        return nu.eye(3)       
        
        
    def __eq__(self,other):
        if isinstance(other,BravaisAtoms):
            # More strict comparison for BravaisAtoms
            # TODO: check additional stuff for other generalized classes
            # (angles, torsions ...)
            if (self.positions-other.positions<1E-13).all() and \
               (self.pbc==other.pbc).all() and \
               self.get_chemical_symbols()==other.get_chemical_symbols() and \
               (self.cell==other.cell).all():
                return True
            else:
                return False
        else:
            # other is probably normal ase.Atoms; more loose check
            if (self.positions-other.positions<1E-13).all() and \
               (self.pbc==other.pbc).all() and \
               self.get_chemical_symbols()==other.get_chemical_symbols() and \
               (self.cell==other.cell).all():
                return True
            else:
                return False
           
        
#    def get_number_of_cells(self):

if __name__=='__main__':
    from ase import *
    atoms=Atoms('C',[(0,0,0)],cell=[(1,1,0),(0,1,0),(0,0,1)])
#    gy=NormalGeometry('C',[(0,0,0)],cell=[(1,1,0),(0,1,0),(0,0,1)])
    gy=NormalGeometry(atoms)
    print gy.vector(array([0,0,0]),n=(1,1,1),r0=[0,0,0],lst=['vec','hat','norm'])
    print gy.vector(array([0,0,0]),n=(1,1,1),lst='dtensor')
    print gy.vector(array([0,0,0]),n=(1,1,1),r0=[2,2,2],lst='vec')
        
        
        
    
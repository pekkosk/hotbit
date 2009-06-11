from ase import Atoms as ase_Atoms
import numpy as nu

class BravaisAtoms(ase_Atoms):
     

    def __init__(self, *args, **kwargs):
        ase_Atoms.__init__(self, *args, **kwargs)
               
    
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
        
    
import numpy as np
from box.mix import phival
from math import sin,cos 
from weakref import proxy

class ContainerTest1:
    
    def __init__(self,atoms,type):
        """
        Container test
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        """
        self.type = 'ContainerTest1'
        assert type==self.type
        self.atoms = proxy(atoms)
        self._set_table() 
        #self.atoms.set_pbc((True,False,False))
        self.atoms.set_pbc((True,True,False))
        
    def __repr__(self):
        return 'NO CONTAINER DOCS'
        
    def set(self,**args):
        #self.atoms.set_pbc((True,False,False))
        self.atoms.set_pbc((True,True,False))
#        if args!={}:
#            print args
#            raise NotImplementedError('Nothing can be set in test container.')

    
    def get_pbc(self):
        """ Return atoms's pbc as normal."""
        #return np.array((True,False,False))
        return np.array((True,True,False))
    
    def get_ase_cell(self):
        """ cell used for visualization """
        return self.atoms.get_cell()
    
    def __eq__(self,other):
        if isinstance(other,ContainerTest1) and np.all(self.get_pbc()==other.get_pbc()) \
           and np.all( (self.atoms.get_cell()-other.atoms.get_cell())<1E-12 ):
            return True
        else:
            return False        
        
        
    def _set_table(self):
        """ Setup group multiplication-like table. """
        eq = (1,0,0)
        self.table = [{'M':np.Inf},{'M':2,'equivalent':eq},{'M':0}]

        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        self._set_table()
        return np.array([[-np.Inf,np.Inf],[-1,1],[0,0]])
    
    
    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        cell = self.atoms.get_cell()
        L = cell[0,0]
        r_mirror = np.array((2*L-r[0],r[1],r[2]))
        dr = r_mirror - r
        rn = r.copy() 
        #if n[1]==-1:
            
        rn = rn + (2*L*n[0],0,0)
        if n[1]==1:
            rn = rn + dr
#            
#        if np.mod(n[0],2)==0:
#            rn[0] += cell[0,0]*n[0]
#        else:
#            rn[0] = (n[0]+1)*cell[0,0] - r[0]  
        return rn
        
    
    def rotation(self,n,angles=False):
        """ No rotation in translations. """
        R = np.eye(3)
        if n[1]==1:
            R[0,0] = -1
        return R
        
#        if np.mod(n[0],2)==0:
#            return np.eye(3)
#        else:
#            R = np.eye(3)
#            R[0,0] = -1
#            return R
    
    


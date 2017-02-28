import numpy as np
from weakref import proxy

class Slab:
    
    def __init__(self,atoms,type):
        """
        Class for Slab lattice, for hotbit.Atoms -class
        
        Symmetry operations are:
        S1) Translation 1 along xy-plane 
        S2) Translation 2 along xy-plane 
        S3) Reflection wrt xy-plane, followed by x1*S1 + x2*S2 translation
            Hence S3**2 = S1**(2*x1)*S2**(2*x2), where x=(x1,x2) and
            xi can be 0 or 0.5.
         
        parameters:
        ===========
        atoms:    hotbit.Atoms -instance
        type:     Should equal to "Slab"           
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        """
        self.type = 'Slab'
        assert type==self.type
        self.atoms = proxy(atoms)
        self.x = np.array([0.,0.])
        atoms.set_pbc(True)
        
    def __repr__(self):
        pbc=self.atoms.get_pbc()
        x=self.x
        cell=self.atoms.get_cell()
        d = []
        for a in range(3):
            d.append( np.linalg.norm(cell[a,:]) )
            
        a12 = np.dot(cell[0],cell[1])/(d[0]*d[1])
        a13 = np.dot(cell[0],cell[2])/(d[0]*d[2])
        a23 = np.dot(cell[1],cell[2])/(d[1]*d[2])
        # the third axis must be orthogonal to two others
        assert np.abs(a13)<1E-6 and np.abs(a23)<1E-6
        x='Slab: pbc:[%i,%i,%i], x=(%.2f,%.1f)' %(pbc[0],pbc[1],pbc[2],x[0],x[1])
        x+='cell:[%.2f,%.2f,%.2f] Ang, cosines(12,13,23):[%.2f,%.2f,%.2f]' %(d[0],d[1],d[2],a12,a13,a23) 
        return x
        
        
    def _set_table(self):
        if np.any( abs(self.x)>1E-12 ):
            x2 = np.array( np.round(2*self.x),int )
            eq = (x2[0],x2[1],0)
            self.table = [{'M':np.Inf},{'M':np.Inf},{'M':2,'equivalent':eq}]
        else:
            self.table = [{'M':np.Inf},{'M':np.Inf},{'M':2}]  
        
        
    def set(self,container=None,x=None):
        """
        Reset x.
         
        parameters:
        ===========
        x:        fractional translation offset related to reflection.
                  Only integers and half-integers allowed.
        """
        if container is not None:
            # copy container
            assert x is None
            self.set(x=container.x)
        else:
            x = np.array(x)
            if x is not None:
                assert np.all( abs(np.round(2*x)-2*x)<1E-15 )
                self.x = np.array(x)
        self._set_table()
        
        
    def get_pbc(self):
        """ Return atoms's pbc as normal."""
        return self.atoms.get_pbc()
    
    def get_ase_cell(self):
        """ cell used for visualization """
        return self.atoms.get_cell()
    
    def __eq__(self,other):
        
        if isinstance(other,Slab) and np.all(self.atoms.get_pbc()==other.atoms.get_pbc()) \
           and np.all(np.abs(self.atoms.get_cell()-other.atoms.get_cell())<1E-12) \
           and np.all(np.abs(self.x-other.x)<1E-12):
            return True
        else:
            return False        
        
        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        return np.array( [[-np.Inf,np.Inf],[-np.Inf,np.Inf],[0,1]] )
        
    
    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        rn=r.copy()
        cell = self.atoms.get_cell()
        for a in range(2):
            rn = rn + n[a]*cell[a,:]
            if n[2]==1:
                rn = rn + cell[a,:]*self.x[a]
        if n[2]==1:
            rn[2] = -rn[2]
        return rn
    
    def rotation(self,n,angles=False):
        """ No rotation in translations. """
        R = np.eye(3)
        if n[2]==0:
            return R
        else:
            R[2,2] = -1.0
            return R
    
    


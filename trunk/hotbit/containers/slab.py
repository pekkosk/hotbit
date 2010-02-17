import numpy as nu
from weakref import proxy

class Slab:
    
    def __init__(self,atoms,type):
        """
        Class for Slab lattice, for hotbit.Atoms -class
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        """
        self.type = 'Slab'
        assert type==self.type
        self.atoms = proxy(atoms)
        pbc = atoms.get_pbc() 
        pbc[2] = True
        atoms.set_pbc(pbc)
        
    def __repr__(self):
        pbc=self.atoms.get_pbc()
        cell=self.atoms.get_cell()
        d = []
        for a in range(3):
            d.append( nu.linalg.norm(cell[a,:]) )
            
            
        a12 = nu.dot(cell[0],cell[1])/(d[0]*d[1])
        a13 = nu.dot(cell[0],cell[2])/(d[0]*d[2])
        a23 = nu.dot(cell[1],cell[2])/(d[1]*d[2])
        # the third axis must be orthogonal to two others
        assert nu.abs(a13)<1E-6 and nu.abs(a23)<1E-6
        x='Slab: pbc:[%i,%i,%i], ' %(pbc[0],pbc[1],pbc[2])
        x+='cell:[%.2f,%.2f,%.2f] Ang, cosines(12,13,23):[%.2f,%.2f,%.2f]' %(d[0],d[1],d[2],a12,a13,a23) 
        return x
        
    def set(self,container=None,**args):
        #print args
        #if args!={}: 
        #    raise NotImplementedError('For Slab use set_pbc and set_cell as normal.')
        pass
    
    def get_pbc(self):
        """ Return atoms's pbc as normal."""
        return self.atoms.get_pbc()
    
    def get_ase_cell(self):
        """ cell used for visualization """
        return self.atoms.get_cell()
    
    def __eq__(self,other):
        if isinstance(other,Slab) and nu.all(self.atoms.get_pbc()==other.atoms.get_pbc()) \
           and nu.all(nu.abs(self.atoms.get_cell()-other.atoms.get_cell())<1E-12):
            return True
        else:
            return False        
        
        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        ranges = []
        pbc = self.atoms.get_pbc()
        for i in range(2):
            if pbc[i]:
                ranges.append([-nu.Inf,nu.Inf])
            else:
                ranges.append([0,0])
        ranges.append([0,1])
        return nu.array(ranges)
    
    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        rn=r.copy()
        cell = self.atoms.get_cell()
        for a in range(2):
            rn = rn + n[a]*cell[a,:]
        if n[2]==1:
            rn[2] = -rn[2]
        return rn
    
    def rotation(self,n):
        """ No rotation in translations. """
        R = nu.eye(3)
        if n[2]==0:
            return R
        else:
            R[2,2] = -1.0
            return R
    
    


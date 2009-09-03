import numpy as nu
from box.mix import phival
from math import sin,cos 
from weakref import proxy

class Bravais:
    
    def __init__(self,atoms,type):
        """
        Class for Bravais lattice, for hotbit.Atoms -class
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        """
        self.type = 'Bravais'
        assert type==self.type
        self.atoms = proxy(atoms)
        
    def __repr__(self):
        pbc=self.atoms.get_pbc()
        cell=self.atoms.get_cell()
        d = []
        for a in range(3):
            d.append( nu.linalg.norm(cell[a,:]) )
            
        a12 = nu.dot(cell[0],cell[1])/(d[0]*d[1])
        a13 = nu.dot(cell[0],cell[2])/(d[0]*d[2])
        a23 = nu.dot(cell[1],cell[2])/(d[1]*d[2])
        x='Bravais: pbc:[%i,%i,%i], ' %(pbc[0],pbc[1],pbc[2])
        x+='cell:[%.2f,%.2f,%.2f] Ang, cosines(12,13,23):[%.2f,%.2f,%.2f]' %(d[0],d[1],d[2],a12,a13,a23) 
        return x
        
    def set(self,**args):
        if args!={}: 
            raise NotImplementedError('For Bravais use set_pbc and set_cell as normal.')
    
    def get_pbc(self):
        """ Return atoms's pbc as normal."""
        return self.atoms.get_pbc()
    
    def get_ase_cell(self):
        """ cell used for visualization """
        return self.atoms.get_cell()
    
    def __eq__(self,other):
        if isinstance(other,Bravais) and all(self.get_pbc()==other.get_pbc()) \
           and self.get_cell()==other.get_cell():
            return True
        else:
            return False        
        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        ranges = []
        pbc = self.atoms.get_pbc()
        for i in range(3):
            if pbc[i]:
                ranges.append([-nu.Inf,nu.Inf])
            else:
                ranges.append([0,0])
        return nu.array(ranges)
    
    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        rn=r.copy()
        cell = self.atoms.get_cell()
        for a in range(3):
            rn = rn + n[a]*cell[a,:]
        return rn
    
    def rotation(self,n):
        """ No rotation in translations. """
        return nu.eye(3)
    
    


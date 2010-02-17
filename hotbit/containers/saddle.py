import numpy as nu
from box.mix import phival
from math import sin,cos 
from weakref import proxy
from box.mix import rotation_matrix
from scipy import linalg

expm = linalg.matfuncs.expm

class Saddle:
    def __init__(self,atoms,type):
        '''
        Class for saddle point structure.
        
        @param: atoms    hotbit.Atoms instance
        @param: type     Should equal to "Saddle"
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        '''
        self.type='Saddle'
        assert type==self.type
        self.atoms = proxy(atoms)
        self.angle1 = None
        self.angle2 = None
        self.n1 = nu.array([0,1,0])
        self.n2 = nu.array([1,0,0])
        self.R = None
        self.M = 50
        
    def __repr__(self):
        x='Saddle: angle1=%.4f, angle2=%.4f, R=%.4f' %(self.angle1,self.angle2,self.R)                                                                                        
        return x
        
            
    def set(self, angle1=None, angle2=None, R=None, container=None):
        """ 
        TODO: doc 
        n1 and n2 should be in xy-plane
        """
        if container!=None:
            assert angle1==None and angle2==None and R==None
            self.set(angle1=container.angle1, angle2=container.angle2, R=container.R)
            
        if angle1!=None:
            self.angle1=angle1
        if angle2!=None:
            self.angle2=angle2
        if R!=None:
            self.R = R
    
        self.atoms.set_pbc((True,True,False))
        self.atoms.set_cell(self.get_ase_cell())  
                   
            
    def __eq__(self,other):
        e=1E-12
        if isinstance(other,Saddle) and \
            abs(self.angle1-other.angle1)<e and \
            abs(self.angle2-other.angle2)<e and \
            abs(self.R-other.R)<e:
            return True
        else:
            return False
    
    def get_ase_cell(self):
        """ cell used for visualization """
        return nu.eye(3)
        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        ranges = nu.array([[-nu.Inf,nu.Inf],[-nu.Inf,nu.Inf],[0,0]])
        return ranges
    
    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        if n[0]==n[1]==0:
            return r.copy()
        A = nu.array( [0,0,2*self.R] )
        
        a1,a2 = self.angle1*n[0], self.angle2*n[1]
        rot1 = rotation_matrix(self.n1,a1/self.M)
        rot2 = rotation_matrix(self.n2,a2/self.M)
        
        rp = r.copy()        
        for i in range(self.M):
            rp = nu.dot(rot1,rp)
            rp = A + nu.dot( rot2,rp-A )
        return rp        
        
    def rotation(self,n):
        """ Rotate around two axes, ordering depending on mode. """
        if n[0]==n[1]==0:
            return nu.eye(3)
        else:
            a1,a2 = self.angle1*n[0], self.angle2*n[1]
            #rot21 = nu.dot( rotation_matrix(self.n2,a2/self.M),rotation_matrix(self.n1,a1/self.M) )
            rot1 = rotation_matrix(self.n1,a1/self.M)
            rot2 = rotation_matrix(self.n2,a2/self.M)
            rot21 = nu.dot(rot2,rot1)
            
            rot = nu.eye(3)       
            for i in range(self.M):
                #rot = nu.dot( rot1,rot )
                #rot = nu.dot( rot2,rot ) 
                rot = nu.dot( rot21,rot )
            return rot
           
        
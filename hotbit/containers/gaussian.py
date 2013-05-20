import numpy as np
from box.mix import phival
from math import sin,cos
from numpy import abs
from weakref import proxy
from box.mix import rotation_matrix

class Gaussian:
    def __init__(self,atoms,type):
        '''
        Class for arbitrary Gaussian curvature.
        
        @param: atoms    hotbit.Atoms instance
        @param: type     Should equal to "Saddle"
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        '''
        self.type='Gaussian'
        assert type==self.type
        self.atoms = proxy(atoms)
        self.angle1 = None
        self.angle2 = None
        self.n1 = np.array([-1,0,0])
        self.n2 = np.array([0,1,0])
        self.R1 = None
        self.R2 = None
        
        raise NotImplementedError('Gaussian container does not work properly')
        
    def __repr__(self):
        x='Gaussian: angle1=%.4f, angle2=%.4f, R1=%.4f, R2=%.4f' %(self.angle1,self.angle2,self.R1,self.R2)                                                                                        
        return x
        
            
    def set(self, angle1=None, angle2=None, R1=None, R2=None, container=None):
        """ 
        TODO: doc 
        n1 and n2 should be in xy-plane
        """
        if container!=None:
            assert angle1==None and angle2==None and R1==None and R2==None
            self.set(angle1=container.angle1, angle2=container.angle2, R1=container.R1, R2=container.R2)
            
        if angle1!=None:
            self.angle1=angle1
        if angle2!=None:
            self.angle2=angle2
        if R1!=None:
            self.R1 = R1
        if R2!=None:
            self.R2 = R2    
    
        self.atoms.set_pbc((True,True,False))
        self.atoms.set_cell(self.get_ase_cell())  
                   
            
    def __eq__(self,other):
        e=1E-12
        if isinstance(other,Saddle) and \
            abs(self.angle1-other.angle1)<e and \
            abs(self.angle2-other.angle2)<e and \
            abs(self.R1-other.R1)<e and \
            abs(self.R2-other.R2)<e:
            return True
        else:
            return False
    
    def get_ase_cell(self):
        """ cell used for visualization """
        return np.eye(3)
        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        ranges = np.array([[-np.Inf,np.Inf],[-np.Inf,np.Inf],[0,0]])
        return ranges
    
    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        if n[0]==n[1]==0:
            return r.copy()
        eps=1E-10
        
        a,b = self.angle1*n[0], self.angle2*n[1]
        R1, R2 = self.R1, self.R2
        R = ( (a*R1)**2+(b*R2)**2 )/( a**2*R1+b**2*R2 )
        
        #print R
        assert abs(R)+eps>=min( abs(R1),abs(R2) )
        A = np.array( [0,0,abs(R1)+abs(R2)] ) # another center of curvature
        rot = self.rotation(n)
        if R>0:
            #return np.array([0,0,0])
            return np.dot(rot,r)
        else:
            if R>0:
                B=(0,0,R1-R)
            else:
                B=A-(0,0,R)
            #B = (0,0,R1+abs(R2)-R)
            #return np.array([0,0,0])
            return B + np.dot(rot,r-B)

    def rotation(self,n,angles=False):
        """ Rotate around two axes, ordering depending on mode. """
        if n[0]==n[1]==0:
            return np.eye(3)
        else:
            R1, R2 = self.R1, self.R2
            a,b = self.angle1*n[0], self.angle2*n[1]
            R = ( (a*R1)**2+(b*R2)**2 )/( a**2*R1+b**2*R2 )

            axis = np.array( [R2*b,R1*a,0] )
            angle = ( a**2*R1+b**2*abs(R2) )/np.sqrt( (a*R1)**2+(b*R2)**2 )
            return rotation_matrix( axis,angle )
        
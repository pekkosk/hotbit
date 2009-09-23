import numpy as nu
from box.mix import phival
from math import sin,cos 
from weakref import proxy
from box.mix import rotation_matrix

class Sphere:
    def __init__(self,atoms,type):
        '''
        Class for spherical boundary conditions.
        
            ______
           /_____/
        
             
        
              |z
              |
              |_________ y
             /
            / 
           /x
        
        @param: atoms    hotbit.Atoms instance
        @param: type     Should equal to "Sphere"
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        '''
        self.type='Sphere'
        assert type==self.type
        self.atoms = proxy(atoms)
        self.angle1 = None
        self.angle2 = None
        self.n1 = None
        self.n2 = None
        self.mode = 4 # mode 4 appears to be the best
        
    def __repr__(self):
        x='Sphere: angle1=%.4f, angle2=%.4f, cos1=(%.2f,%.2f,%.2f), cos2=(%.2f,%.2f,%.2f)' %(self.angle1,self.angle2,self.n1[0],self.n1[1],self.n1[2],self.n2[0],self.n2[1],self.n2[2])                                                                                        
        return x
        
            
    def set(self, angle1=None, angle2=None, n1=None, n2=None, mode=None, container=None):
        """ 
        TODO: doc 
        n1 and n2 should be in xy-plane
        """
        if container!=None:
            assert angle1==None and angle2==None and n1==None and n2==None and mode==None
            self.set(angle1=container.angle1, angle2=container.angle2,
                     n1=container.n1, n2=container.n2, mode=container.mode)
            
        if angle1!=None:
            self.angle1=angle1
        if angle2!=None:
            self.angle2=angle2
        if n1!=None:
            assert abs(n1[2])<1E-10
            if self.n1!=None: 
                raise AssertionError('Rotation axis n1 cannot be changed.')
            self.n1 = n1
        if n2!=None:
            assert abs(n2[2])<1E-10
            if self.n2!=None: 
                raise AssertionError('Rotation axis n2 cannot be changed.')
            self.n2 = n2
        if mode!=None:
            self.mode = mode
    
        if self.n1!=None:
            self.n1 = self.n1/nu.linalg.norm(self.n1)
        if self.n2!=None:
            self.n2 = self.n2/nu.linalg.norm(self.n2)
        self.atoms.set_pbc((True,True,False))
        self.atoms.set_cell(self.get_ase_cell())  
                   
            
    def __eq__(self,other):
        e=1E-12
        if isinstance(other,Sphere) and \
            abs(self.angle1-other.angle1)<e and \
            abs(self.angle2-other.angle2)<e and \
            nu.linalg.norm(self.n1-other.n1)<e and \
            nu.linalg.norm(self.n2-other.n2)<e:
            return True
        else:
            return False
    
    def get_ase_cell(self):
        """ cell used for visualization """
        h = max(self.atoms.get_positions()[:,2])
        phi1 = phival(self.n1[0],self.n1[1])
        phi2 = phival(self.n2[0],self.n2[1])
        a1,l1 = phi1-nu.pi/2, h*self.angle1
        a2,l2 = phi2-nu.pi/2, h*self.angle2
        return nu.array( [[l1*cos(a1),l1*sin(a1),0],[l2*cos(a2),l2*sin(a2),0],[0,0,h]] )
        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        ranges = nu.array([[-nu.Inf,nu.Inf],[-nu.Inf,nu.Inf],[0,0]])
        return ranges
    
    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        return nu.dot( self.rotation(n),r )

    def rotation(self,n):
        """ Rotate around two axes, ordering depending on mode. """
        if n[0]==n[1]==0:
            return nu.eye(3)
        if self.mode==1 or self.mode==2:
            R1n = rotation_matrix( self.n1,n[0]*self.angle1 )
            R2n = rotation_matrix( self.n2,n[1]*self.angle2 )
            if self.mode==1:
                return nu.dot(R2n,R1n)
            else:
                return nu.dot(R1n,R2n)
        elif self.mode==3:
            # rotate R2**l2 * R1**l1 * (R2*R2)**m
            na = (abs(n[0]),abs(n[1]),0)
            m = min(na[:2])            
            R1 = rotation_matrix( self.n1,nu.sign(n[0])*self.angle1 )
            R2 = rotation_matrix( self.n2,nu.sign(n[1])*self.angle2 )
            R21 = nu.dot(R2,R1)
            R = nu.eye(3)
            for i in range(m):
                R = nu.dot(R,R21)
            # now rotate with the 'left-over'
            l1 = na[0]-m
            l2 = na[1]-m
            for i in range(l1):
                R = nu.dot(R1,R)
            for i in range(l2):
                R = nu.dot(R2,R)
            return R
        elif self.mode==4:
            axis = n[0]*self.n1 + n[1]*self.n2
            a1,a2 = self.angle1*n[0], self.angle2*n[1]
            angle = nu.sqrt( a1**2 + a2**2 + 2*a1*a2*nu.dot(self.n1,self.n2) )
            return rotation_matrix( axis,angle )
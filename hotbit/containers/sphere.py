import numpy as np
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
        self.par = {'angle1':(0,1),'angle2':(0,2),'n1':(1,slice(0,2)),'n2':(2,slice(0,2)),'mode':(1,2)}
        self._set_table()
        
    def _set_table(self): 
        self.table = [{'M':np.Inf},{'M':np.Inf},{'M':1}]
     
    def get_type(self):
        return self.type
        
    def __repr__(self):
        a1,a2,n1,n2 = [self.get(k) for k in ['angle1','angle2','n1','n2']]
        x='Sphere: angle1=%.4f, angle2=%.4f, cos1=(%.2f,%.2f), cos2=(%.2f,%.2f)' %(a1,a2,n1[0],n1[1],n2[0],n2[1])                                                                                        
        return x
    
    def get_table(self):
        return [{'M':np.Inf},{'M':np.Inf},{'M':1}]
    
    def get(self,key):
        """
        Get container parameters.
        key = 
        """
        cell = self.atoms.get_cell()
        x = cell[self.par[key]]
        if key=='mode':
            return int(x)  
        elif key in ['n1','n2']:
            return np.concatenate((x,[0]))
        else:
            return x    
    
    def _set(self,**kwargs):
        """
        Save n1 and n2 only as two-dimensional vectors (ignore z-component)
        """
        assert len(kwargs)==1
        for key in kwargs:
            cell = self.atoms.get_cell()
            value = kwargs[key]
            if key in ['n1', 'n2']:
                value = value[0:2]
            cell[self.par[key]] = value
            self.atoms.set_cell(cell)
            
    def set(self, angle1=None, angle2=None, n1=None, n2=None, mode=None, container=None):
        """ 
        Rotate angle1 around axis n1, rotate angle2 around axis n1 (wrt origin).
        n1 and n2 should be in xy-plane
        """
        if container!=None:
            assert angle1==None and angle2==None and n1==None and n2==None and mode==None
            self.set(angle1=container.angle1, angle2=container.angle2,
                     n1=container.n1, n2=container.n2, mode=container.mode)
            
        if angle1!=None:
            self._set(angle1=angle1)
        if angle2!=None:
            self._set(angle2=angle2)
        if n1!=None:
            assert abs(n1[2])<1E-10
            if self.n1!=None: 
                raise AssertionError('Rotation axis n1 cannot be changed.')
            self._set(n1 = n1)
        if n2!=None:
            assert abs(n2[2])<1E-10
            if self.n2!=None: 
                raise AssertionError('Rotation axis n2 cannot be changed.')
            self._set(n2 = n2)
        if mode!=None:
            self._set(mode = mode)
    
        if self.n1!=None:
            self._set(n1 = self.n1/np.linalg.norm(self.n1))
        if self.n2!=None:
            self._seT(n2 = self.n2/np.linalg.norm(self.n2))
        self.atoms.set_pbc((True,True,False))
        #self.atoms.set_cell(self.get_ase_cell())  
                   
            
    def __eq__(self,other):
        return self.atoms == other.atoms
    #===========================================================================
    # def __eq__(self,other):
    #    e=1E-12
    #    if isinstance(other,Sphere) and \
    #        abs(self.angle1-other.angle1)<e and \
    #        abs(self.angle2-other.angle2)<e and \
    #        np.linalg.norm(self.n1-other.n1)<e and \
    #        np.linalg.norm(self.n2-other.n2)<e:
    #        return True
    #    else:
    #        return False
    # 
    #===========================================================================
    def get_ase_cell(self):
        """ cell used for visualization """
        h = max(self.atoms.get_positions()[:,2])
        phi1 = phival(self.n1[0],self.n1[1])
        phi2 = phival(self.n2[0],self.n2[1])
        a1,l1 = phi1-np.pi/2, h*self.angle1
        a2,l2 = phi2-np.pi/2, h*self.angle2
        return np.array( [[l1*cos(a1),l1*sin(a1),0],[l2*cos(a2),l2*sin(a2),0],[0,0,h]] )
        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        ranges = np.array([[-np.Inf,np.Inf],[-np.Inf,np.Inf],[0,0]])
        return ranges
    
    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        return np.dot( self.rotation(n),r )

    def rotation(self,n,angles=False):
        """ Rotate around two axes, ordering depending on mode. """
        if n[0]==n[1]==0:
            return np.eye(3)
        angle1, angle2, n1, n2, mode = [self.get(k) for k in ['angle1','angle2','n1','n2','mode']]
        if angles and mode!=4:
            raise NotImplementedError('Sphere rotation implemented only for mode 4 so far.')
        if self.mode==1 or self.mode==2:
            R1n = rotation_matrix( n1,n[0]*angle1 )
            R2n = rotation_matrix( n2,n[1]*angle2 )
            if mode==1:
                return np.dot(R2n,R1n)
            else:
                return np.dot(R1n,R2n)
        elif mode==3:
            # rotate R2**l2 * R1**l1 * (R2*R2)**m
            na = (abs(n[0]),abs(n[1]),0)
            m = min(na[:2])            
            R1 = rotation_matrix( n1,np.sign(n[0])*angle1 )
            R2 = rotation_matrix( n2,np.sign(n[1])*angle2 )
            R21 = np.dot(R2,R1)
            R = np.eye(3)
            for i in range(m):
                R = np.dot(R,R21)
            # now rotate with the 'left-over'
            l1 = na[0]-m
            l2 = na[1]-m
            for i in range(l1):
                R = np.dot(R1,R)
            for i in range(l2):
                R = np.dot(R2,R)
            return R
        elif self.mode==4:
            axis = angle1*n[0]*n1 + angle2*n[1]*n2
            a1,a2 = angle1*n[0], angle2*n[1]
            angle = np.sqrt( a1**2 + a2**2 + 2*a1*a2*np.dot(n1,n2) )
            if angles:
                axis = axis/np.linalg.norm(axis)
                return np.pi/2.,np.arctan2(axis[1],axis[0]),angle
            else:
                return rotation_matrix( axis,angle )
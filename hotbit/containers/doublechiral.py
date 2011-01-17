import numpy as np
from box.mix import phival
from math import sin,cos 
from weakref import proxy


class DoubleChiral:
    
    def __init__(self,atoms,type):
        '''
        Class for doubly chiral boundary conditions.
        
        Symmetry operations are:
        1) 180 degrees rotation around z-axis + 
           completing chiral transformation if system has translational symmetry
           
           S_1: R_z(n_1*pi + n_1*angle*x) + n_1*x*height*z
           
        3) Rotation around z-axis and simultaneous translation along z-axis:
        
           S_3: R_z(n_3*angle)r + n_3*height*z
           
        The total symmetry operation is hence
        
        S(n_1,n_3): r' = R_z(n_3*angle + n_1*(pi + angle*x) + 
                         (n_3 + n_1*x)*height*z 
                         
                    where n_3 = 0,+-1,+-2,+-3,...
                          n_1 = 0,1
                            x = 0....1 (fraction of full translation)
        
        parameters:
        ===========
        
        atoms:    hotbit.Atoms -instance
        type:     Should equal to "DoubleChiral"
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        '''
        self.type='DoubleChiral'
        assert type==self.type
        self.atoms = proxy(atoms)
        self.par = {'angle':(1,0),'height':(2,0),'x':(0,1)}
        self.atoms.set_pbc((True,False,True))
        
        
    def __repr__(self):
        height, angle, x = self.get('height'), self.get('angle'), self.get('x')
        if angle<1E-14:
            st = 'DoubleChiral: angle=0.0, height=%.4f Ang, x=%.4f' %(height,x)
        else:
            st = 'DoubleChiral: angle=%.4f (2*pi/%.2f), height=%.4f Ang, x=%.4f' %(angle,2*np.pi/angle,height,x)
        return st
    
    
    def get_table(self):
        if abs(self.get('x'))>1E-12:
            eq = ( 0,0,int(np.round(2*self.get('x'))) )
            return [{'M':2,'equivalent':eq},{'M':1},{'M':np.Inf}]
        else:
            return [{'M':2},{'M':1},{'M':np.Inf}]
    
    
    def get(self,key):
        """
        Return current angle, height, or x
        """
        cell = self.atoms.get_cell()
        if key in ['angle','height','x']:
            return cell[self.par[key]]
        else:
            raise AssertionError('Invalid keyword %s' %key)
    

    def _set(self,**kwargs):
        assert len(kwargs)==1
        for key in kwargs:
            cell = self.atoms.get_cell()
            cell[self.par[key]] = kwargs[key]
            self.atoms.set_cell(cell)

        
    def set(self,angle=None,height=None,x=None,scale_atoms=False,container=None):
        """
        Reset angle, height, or x, and maybe scale atoms.
         
        parameters:
        ===========
        height:   Height of the primitive cell in z-direction
        angle:    angle (in radians) of rotation
        x:        fractional translation offset related to 180 rotation
                  Only integers and half-integers allowed.
        """
        if container!=None:
            # copy container
            assert angle==None and height==None and x==None and scale_atoms==False
            self.set(angle=container.get('angle'),height=container.get('height'),x=container.get('x'))
        else:
            if x!=None:
                assert abs(np.round(2*x)-2*x)<1E-15
            if not scale_atoms:
                if angle!=None: self._set(angle=angle)
                if height!=None: self._set(height=height)
                if x!=None: self._set(x=x)
            else:
                if x!=None:
                    raise AssertionError('It is probably illegal to change x. This changes the symmetry, right?')
                if angle is None:
                    da = 0.0
                else:
                    da = angle - self.get('angle')
                    self._set(angle=angle)
                    
                old_height = self.get('height')
                if height != None:
                    self._set(height=height)
                    
                newr = []
                for r in self.atoms.get_positions():
                    x,y = r[0],r[1]
                    rad = np.sqrt( x**2+y**2 )
                    frac = r[2]/old_height
                    # twist atoms z/h * da (more)
                    newphi = phival(x,y) + frac * da
                    newr.append( [rad*np.cos(newphi),rad*np.sin(newphi),frac*self.get('height')] )
                self.atoms.set_positions(newr)     
        
    def __eq__(self,other):
        return self.atoms == other.atoms
            
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        return np.array( [[0,1],[0,0],[-np.Inf,np.Inf]] )
    
    def transform(self,r,n):
        """ See init doc for transformation."""
        R = self.rotation(n)
        return np.dot(R,r) + (0,0,(n[2]+n[0]*self.get('x'))*self.get('height'))
                         
                         
    def rotation(self,n):
        """ Return the (active) rotation matrix for symmetry operation n. """
        angle = n[2]*self.get('angle') + n[0]*(np.pi+self.get('angle')*self.get('x'))
        R = np.array([[cos(angle),-sin(angle),0],[sin(angle),cos(angle),0],[0,0,1]])
        return R 
    

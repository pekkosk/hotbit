import numpy as np
from box.mix import phival
from math import sin,cos 
from weakref import proxy


class Chiral:
    
    def __init__(self,atoms,type):
        '''
        Class for chiral boundary conditions.
        
        @param: atoms    hotbit.Atoms -instance
        @param: type     Should equal to "Chiral"
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        '''
        self.type='Chiral'
        assert type == self.type
        self.atoms = proxy(atoms)
        self.par = {'height':(1,0),'angle':(2,0)}
        self.atoms.set_pbc( (False,False,True) )
        
    def get_type(self):
        return self.type
        
    def get_table(self): 
        return [{'M':1},{'M':1},{'M':np.Inf}]
        
    def __repr__(self):
        angle = self.get('angle')
        height = self.get('height')
        if angle==None:
            raise AssertionError('Chiral angle was not set yet.')
        if angle<1E-14:
            x='Chiral: angle=0.0, height=%.4f Ang' %(height)
        else:
            x='Chiral: angle=%.4f (2*pi/%.2f), height=%.4f Ang' %(angle,2*np.pi/angle,height)
        return x
    
    
    def get(self,key):
        """
        Return container parameters.
        
        parameters:
        ===========
        key:    'angle','height'
        """
        if key not in ['angle','height']:
            raise AssertionError('Invalid keyword %s' %key)
        if key=='angle':
            return self.atoms.get_cell()[self.par['angle']]
        elif key=='height':
            return self.atoms.get_cell()[self.par['height']]             
 
    def _set(self,**kwargs):
        for key in kwargs:
            cell = self.atoms.get_cell()
            cell[self.par[key]] = kwargs[key]
            self.atoms.set_cell(cell)
 
        
    def set(self,angle=None,height=None,scale_atoms=False,container=None):
        """
        Reset angle or height, and maybe scale atoms.
         
        @param: height   Height of the primitive cell in z-direction
        @param: angle    angle (in radians) of rotation
        """
        if container!=None:
            # copy container
            assert angle==None and height==None and scale_atoms==False
            self.set( angle=container.get('angle'),height=container.get('height') )
        else:
            if not scale_atoms:
                if angle!=None: self._set(angle=angle)
                if height!=None: self._set(height=height)
            else:
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
        return np.array( [[0,0],[0,0],[-np.Inf,np.Inf]] )
    
    def transform(self,r,n):
        """ Rotate r by n2*angle + translate (in z) by n2*height."""
        R = self.rotation(n)
        return np.dot(R,r) + (0,0,n[2]*self.get('height'))
    
    def rotation(self,n):
        """ Return the (active) rotation matrix for symmetry operation n. """
        angle = n[2]*self.get('angle')
        R = np.array([[cos(angle),-sin(angle),0],[sin(angle),cos(angle),0],[0,0,1]])
        return R 
    

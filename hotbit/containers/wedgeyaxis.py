import numpy as np
from box.mix import phival
from math import sin,cos 
from weakref import proxy
import warnings

class WedgeYAxis:
    
    def __init__(self,atoms,type):
        '''
        Class for wedge boundary conditions.
        
        Rotation is around *** y-axis ***
        
           ^ x-axis
           |    /
           |   /
           |  /  
           | /  angle 
           +----------> z-axis
        
        @param: atoms    hotbit.Atoms instance
        @param: type     Should equal to "CustomContainer"
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        '''
        self.type='WedgeYAxis'
        assert type==self.type
        self.atoms = proxy(atoms)
        self.par = {'height':(1,0),'angle':(2,0),'pbc':(0,1),'physical':(1,1)}
        self.atoms.set_pbc((True,False,atoms.get_pbc()[2]))

      
    def get_type(self):
        return self.type  
        
    def __repr__(self):
        angle, height, pbc, physical = self.get('angle'), self.get('height'), self.get('pbc'), self.get('physical')
        x='Wedge: angle=%.4f (2*pi/%.2f, ' %(angle,2*np.pi/angle)
        if physical:
            x+='physical), '
        else:
            x+='not physical), '
        x+='height=%.4f Ang ' %height
        if pbc:
            x+='(y:pbc)'
        else:
            x+='(y:no pbc)'
        return x
    
    def get_table(self):
        M = int( round(2*np.pi/self.get('angle')) )
        if self.get('pbc'):
            return [{'M':M},{'M':np.Inf},{'M':1}]
        else:
            return [{'M':M},{'M':1},{'M':1}]
        
            
    def get(self,key):
        """
        Get container parameters
        
        key: 'angle','height','pbc','physical'
        """
        if key=='pbc':
            return self.atoms.get_pbc()[1]
        else:
            x = self.atoms.get_cell()[self.par[key]]
            if key in ['angle','height']:
                return x  
            else:
                return bool(np.round(x)) 
        
    def _set(self,**kwargs):
        assert len(kwargs)==1
        if 'pbc' in kwargs:
            self.atoms.set_pbc( (True,kwargs['pbc'],False) )
        else:
            for key in kwargs:
                cell = self.atoms.get_cell()
                cell[self.par[key]] = kwargs[key]
                self.atoms.set_cell(cell)
                
            
    def set(self, angle=None, height=None, M=None, physical=True, pbc=None, scale_atoms=False, container=None):
        """ Only height can be reset, not angle. 

        parameters:
        ===========
        angle    angle (in radians) of the wedge (and M=None)
        height   Height of the primitive cell in z-direction
        M        set angle to 2*pi/M (and angle=None)
        physical (only if M=None) if angle is small, it does not be
                         exactly 2*pi/integer, i.e. situation has no physical meaning
                         (use for calculating stuff continuously)
        pbc     True if wedge is periodic in y-direction
        scale_atoms Scale atoms according to changes in parameters 
        """
        if container!=None:
            assert angle==None and height==None and M==None and pbc==None
            self.set(angle=container.get('angle'),height=container.get('height'),\
                     physical=container.get('physical'), pbc=container.atoms.get_pbc()[2])
        
        if angle!=None or M!=None:
            #assert not scale_atoms
            assert not (angle!=None and M!=None)
            old_angle = self.get('angle')
            if M != None:
                assert isinstance(M,int)
                self._set(angle=2*np.pi/M)
            elif angle != None:
                M = int( round(2*np.pi/angle) )
                self._set(angle=angle)
                
            # check parameters
            self._set( physical=float(physical) )
            if self.get('angle')<1E-6:
                raise Warning('Too small angle (%f) may bring numerical problems.' %self.get('angle'))
            if self.get('angle')>np.pi:
                raise AssertionError('angle>pi')
            if np.abs(M-2*np.pi/self.get('angle'))>1E-12 and self.get('physical'): 
                raise AssertionError('angle not physical: angle != 2*pi/M')
            if not self.get('physical') and M<20:
                warnings.warn('Quite large, non-physical angle 2*pi/%.4f.' %(2*np.pi/self.get('angle')) )
            
            if scale_atoms:
                if abs(old_angle)<1E-10:
                    raise ValueError('Atoms cannot be scaled; old wedge angle too small.')
                newr = []
                for r in self.atoms.get_positions():
                    z,x = r[2],r[0]
                    rad = np.sqrt( z**2+x**2 )
                    newphi = phival(z,x)*(self.get('angle')/old_angle)
                    newr.append( [rad*np.sin(newphi),r[1],rad*np.cos(newphi)] )
                self.atoms.set_positions(newr)
                
        if height!=None:
            if scale_atoms:
                r = self.atoms.get_positions()
                r[:,1] = r[:,1] * height/self.get('height')
                self.atoms.set_positions(r)
            self._set(height=height)
            
        if pbc!=None:
            self._set(pbc=float(pbc))
              
            
    def __eq__(self,other):
        return self.atoms == other.atoms
        
        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        M = int( round(2*np.pi/self.get('angle')) )
        i = M/2
        zi = 0
        if self.get('pbc'):
            zi = np.Inf
        if np.mod(M,2)==1:
            ranges = np.array([[-i,i],[-zi,zi],[0,0]])
        else:
            ranges = np.array([[-i+1,i],[-zi,zi],[0,0]])
        return ranges
    
    
    def transform(self,r,n):
        """ Rotate r by n2*angle. """
        R = self.rotation(n)
        trans = np.zeros((3))
        if self.get('pbc'):
            trans = n[1]*np.array([0,self.get('height'),0])
        return np.dot(R,r) + np.array(trans)
    
    def rotation(self,n,angles=False):
        """ Active rotation matrix of given angle wrt. y-axis."""
        angle = n[0]*self.get('angle')
        if angles:
            return np.pi/2,np.pi/2,angle 
        else:
            #R = np.array([[cos(angle),-sin(angle),0],[sin(angle),cos(angle),0],[0,0,1]])
            R = np.array([[cos(angle),0,sin(angle)],[0,1,0],[-sin(angle),0,cos(angle)]])
            return R


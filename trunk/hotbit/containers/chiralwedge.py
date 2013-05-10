import numpy as np
from box.mix import phival
from math import sin,cos 
from weakref import proxy
import warnings

class ChiralWedge:
    
    def __init__(self,atoms,type):
        '''
        Class for chiral+wedge boundary conditions.
               
        @param: atoms    hotbit.Atoms instance
        @param: type     Should equal to "ChiralWedge"
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        '''
        self.type='ChiralWedge'
        assert type==self.type
        self.atoms = proxy(atoms)
        self.par = {'height':(1,0),'twist':(0,1),'angle':(2,0),'physical':(1,1)}        
        self.atoms.set_pbc((True,False,True))
        #self._set_table()
      
    def get_type(self):
        return self.type  
        
    def __repr__(self):
        twist, angle, height, physical = self.get('twist'), self.get('angle'), self.get('height'), self.get('physical')
        x='ChiralWedge: angle=%.4f (2*pi/%.2f, ' %(angle,2*np.pi/angle)
        if physical:
            x+='physical), '
        else:
            x+='not physical), '
        x+='height=%.4f Ang ' %height
        x+='twist angle %.4f' %twist
        return x
    
    def get_table(self):
        M = int( round(2*np.pi/self.get('angle')) )
        return [{'M':M},{'M':1},{'M':np.Inf}]
                    
    def get(self,key):
        """
        Get container parameters
        
        key: 'angle','height','twist','physical'
        """
        x = self.atoms.get_cell()[self.par[key]]
        if key in ['angle','height','twist']:
            return x  
        else:
            return bool(np.round(x)) 
        
    def _set(self,**kwargs):
        assert len(kwargs)==1
        for key in kwargs:
            cell = self.atoms.get_cell()
            cell[self.par[key]] = kwargs[key]
            self.atoms.set_cell(cell)
                
            
    def set(self, angle=None, height=None, M=None, physical=True, twist=None, scale_atoms=False, container=None):
        """ 
        
        parameters:
        ===========
        angle    angle (in radians) of the wedge (and M=None)
        height   Height of the primitive cell in z-direction
        M        set angle to 2*pi/M (and angle=None)
        physical (only if M=None) if angle is small, it does not be
                         exactly 2*pi/integer, i.e. situation has no physical meaning
                         (use for calculating stuff continuously)
        twist     The twist angle for z-translation
        scale_atoms Scale atoms according to changes in parameters 
        """
        if container!=None:
            assert angle==None and height==None and M==None and twist==None
            self.set(angle=container.get('angle'),height=container.get('height'),\
                     physical=container.get('physical'), twist=container.get('twist'))
        
        if angle!=None or M!=None:
            #assert not scale_atoms
            assert not (angle!=None and M!=None)
            old_angle = self.get('angle')
            if M != None:
                assert isinstance(M,int)
                self._set(angle=2*np.pi/M)
            elif angle != None:
                M = np.abs(int( round(2*np.pi/angle) ))
                self._set(angle=angle)
                
            # check parameters
            self._set( physical=float(physical) )
            if np.abs(self.get('angle'))<1E-6:
                raise Warning('Too small angle (%f) may bring numerical problems.' %self.get('angle'))
            if self.get('angle')>np.pi:
                raise AssertionError('angle>pi')
            if np.abs(M-2*np.pi/np.abs(self.get('angle')))>1E-12 and self.get('physical'): 
                raise AssertionError('angle not physical: angle != 2*pi/M')
            if not self.get('physical') and M<20:
                warnings.warn('Quite large, non-physical angle 2*pi/%.4f.' %(2*np.pi/self.get('angle')) )
            
            if scale_atoms:
                if abs(old_angle)<1E-10:
                    raise ValueError('Atoms cannot be scaled; old wedge angle too small.')
                newr = []
                for r in self.atoms.get_positions():
                    x,y = r[0],r[1]
                    rad = np.sqrt( x**2+y**2 )
                    newphi = phival(x,y)*(self.get('angle')/old_angle)
                    newr.append( [rad*np.cos(newphi),rad*np.sin(newphi),r[2]] )
                self.atoms.set_positions(newr)
                
        if height!=None:
            if scale_atoms:
                r = self.atoms.get_positions()
                r[:,2] = r[:,2] * height/self.get('height')
                self.atoms.set_positions(r)
            self._set(height=height)
            
        if twist!=None:
            if scale_atoms:
                raise NotImplementedError('Atom rescale with twist not implemented.')
            self._set(twist=twist)
            
        #self._set_table()  
              
            
    def __eq__(self,other):
        return self.atoms == other.atoms
        
        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        M = int( round(2*np.pi/np.abs(self.get('angle'))) )
        i = M/2
        zi = 0
        if np.mod(M,2)==1:
            ranges = np.array([[-i,i],[0,0],[-np.Inf,np.Inf]])
        else:
            ranges = np.array([[-i+1,i],[0,0],[-np.Inf,np.Inf]])
        return ranges
    
    
    def transform(self,r,n):
        """ Rotate around z r by (n2*angle+n0*twist) and translate by n0*height. """
        R = self.rotation(n)
        trans = np.zeros((3))
        trans = n[2]*np.array([0,0,self.get('height')])
        return np.dot(R,r) + np.array(trans)
    
    def rotation(self,n):
        """ Active rotation matrix of given angle wrt. z-axis."""
        angle = n[0]*self.get('angle') + n[2]*self.get('twist')
        R = np.array([[cos(angle),-sin(angle),0],[sin(angle),cos(angle),0],[0,0,1]])
        return R


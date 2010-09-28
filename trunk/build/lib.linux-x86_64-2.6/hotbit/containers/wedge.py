import numpy as np
from box.mix import phival
from math import sin,cos 
from weakref import proxy
import warnings

class Wedge:
    
    def __init__(self,atoms,type):
        '''
        Class for wedge boundary conditions.
        
           ^ y-axis
           |    /
           |   /
           |  /  
           | /  angle 
           +----------> x-axis
        
        @param: atoms    hotbit.Atoms instance
        @param: type     Should equal to "Wedge"
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        '''
        self.type='Wedge'
        assert type==self.type
        self.atoms = proxy(atoms)
        self.height = None
        self.angle = None
        self.physical = None
        self.pbcz = False
        
        
    def __repr__(self):
        x='Wedge: angle=%.4f (2*pi/%.2f, ' %(self.angle,2*np.pi/self.angle)
        if self.physical:
            x+='physical), '
        else:
            x+='not physical), '
        x+='height=%.4f Ang ' %self.height
        if self.atoms.get_pbc()[2]:
            x+='(z:pbc)'
        else:
            x+='(z:no pbc)'
        return x
    
    
    def _set_table(self):
        if self.pbcz:
             self.table = [{'M':self.M},{'M':1},{'M':np.Inf}]
        else:
             self.table = [{'M':self.M},{'M':1},{'M':1}]
            
            
    def set(self, angle=None, height=None, M=None, physical=True, pbcz=None, scale_atoms=False, container=None):
        """ Only height can be reset, not angle. 

        @param: angle    angle (in radians) of the wedge (and M=None)
        @param: height   Height of the primitive cell in z-direction
        @param: M        set angle to 2*pi/M (and angle=None)
        @param: physical (only if M=None) if angle is small, it does not be
                         exactly 2*pi/integer, i.e. situation has no physical meaning
                         (use for calculating stuff continuously)
        @param: pbcz     True if wedge is periodic in z-direction
        @param: scale_atoms Scale atoms according to changes in parameters 
        """
        if container!=None:
            assert angle==None and height==None and M==None and pbcz==None
            self.set(angle=container.angle,height=container.height,\
                     physical=container.physical, pbcz=container.atoms.get_pbc()[2])
        
        if angle!=None or M!=None:
            #assert not scale_atoms
            assert not (angle!=None and M!=None)
            if self.angle==None and scale_atoms:
                raise AssertionError('Atoms cannot be scaled; angle was not set yet.')
            old_angle = self.angle
            if M != None:
                assert isinstance(M,int)
                self.angle = 2*np.pi/M
                self.M = M
            elif angle != None:
                self.M = int( round(2*np.pi/angle) )
                self.angle = angle
                
            # check parameters
            self.physical = physical
            if self.angle<1E-6:
                raise Warning('Too small angle (%f) may bring numerical problems.' %self.angle)
            if self.angle>np.pi:
                raise AssertionError('angle>pi')
            if np.abs(self.M-2*np.pi/self.angle)>1E-12 and physical: 
                raise AssertionError('angle not physical: angle != 2*pi/M')
            if not physical and self.M<20:
               warnings.warn('Quite large, non-physical angle 2*pi/%.4f.' %(2*np.pi/self.angle) )
            
            if scale_atoms:
                if abs(old_angle)<1E-10:
                    raise ValueError('Atoms cannot be scaled; old wedge angle too small.')
                newr = []
                for r in self.atoms.get_positions():
                    x,y = r[0],r[1]
                    rad = np.sqrt( x**2+y**2 )
                    newphi = phival(x,y)*(self.angle/old_angle)
                    newr.append( [rad*np.cos(newphi),rad*np.sin(newphi),r[2]] )
                self.atoms.set_positions(newr)
                
        if height!=None:
            if scale_atoms:
                r = self.atoms.get_positions()
                r[:,2] = r[:,2] * height/self.height
                self.atoms.set_positions(r)
            self.height = height
            self.atoms.set_cell( self.get_ase_cell() )
            
        if pbcz!=None:
            self.pbcz = pbcz
            
        self.atoms.set_pbc((True,False,self.pbcz))
        self.atoms.set_cell(self.get_ase_cell())
        self._set_table()  
              
            
    def __eq__(self,other):
        if isinstance(other,Wedge) and abs(self.height-other.height)<1E-12 \
           and abs(self.angle-other.angle)<1E-12 and np.all(self.atoms.get_pbc()==other.atoms.get_pbc())\
           and self.physical==other.physical:
            return True
        else:
            return False
    
    def get_ase_cell(self):
        """ cell used for visualization """
        l = max(self.atoms.get_positions()[:,0])*1.5
        return np.array( [[l,0,0],[l*cos(self.angle),l*sin(self.angle),0],[0,0,self.height]] )
        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        if self.height==None or self.angle==None or self.physical==None:
            raise RuntimeError("Wedge's angle or height is not set yet.")
        i = self.M/2
        zi = 0
        if self.pbcz:
            zi = np.Inf
        if np.mod(self.M,2)==1:
            ranges = np.array([[-i,i],[0,0],[-zi,zi]])
        else:
            ranges = np.array([[-i+1,i],[0,0],[-zi,zi]])
        return ranges
    
    def transform(self,r,n):
        """ Rotate r by n2*angle. """
        R = self.rotation(n)
        trans = np.zeros((3))
        if self.pbcz:
            trans = n[2]*np.array([0,0,self.height])
        return np.dot(R,r) + np.array(trans)
    
    def rotation(self,n):
        """ Active rotation matrix of given angle wrt. z-axis."""
        angle = n[0]*self.angle
        R = np.array([[cos(angle),-sin(angle),0],[sin(angle),cos(angle),0],[0,0,1]])
        return R


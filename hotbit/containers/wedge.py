import numpy as nu
from box.mix import phival
from math import sin,cos 
from weakref import proxy

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
        x='Wedge: angle=%.4f (2*pi/%.2f, ' %(self.angle,2*nu.pi/self.angle)
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
                self.angle = 2*nu.pi/M
                self.M = M
            elif angle != None:
                self.M = int( round(2*nu.pi/angle) )
                self.angle = angle
                
            # check parameters
            self.physical = physical
            if self.angle<1E-6:
                raise Warning('Too small angle (%f) may bring numerical problems.' %self.angle)
            if self.angle>nu.pi:
                raise AssertionError('angle>pi')
            if nu.abs(self.M-2*nu.pi/self.angle)>1E-12 and physical: 
                raise AssertionError('angle not physical: angle != 2*pi/M')
            if not physical and self.M<20:
                raise AssertionError('Too large, non-physical angle.')
            
            if scale_atoms:
                if abs(old_angle)<1E-10:
                    raise ValueError('Atoms cannot be scaled; old wedge angle too small.')
                newr = []
                for r in self.atoms.get_positions():
                    x,y = r[0],r[1]
                    rad = nu.sqrt( x**2+y**2 )
                    newphi = mix.phival(x,y)*(self.angle/old_angle)
                    newr.append( [rad*nu.cos(newphi),rad*nu.sin(newphi),r[2]] )
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
                   
            
    def __eq__(self,other):
        if isinstance(other,Wedge) and abs(self.height-other.height)<1E-12 \
           and abs(self.angle-other.angle)<1E-12 and all(self.atoms.get_pbc()==other.atoms.get_pbc())\
           and self.physical==other.physical:
            return True
        else:
            return False
    
    def get_ase_cell(self):
        """ cell used for visualization """
        l = max(self.atoms.get_positions()[:,0])*1.5
        return nu.array( [[l,0,0],[l*cos(self.angle),l*sin(self.angle),0],[0,0,self.height]] )
        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        if self.height==None or self.angle==None or self.physical==None:
            raise RuntimeError("Wedge's angle or height is not set yet.")
        i = self.M/2
        if nu.mod(self.M,2)==1:
            ranges = nu.array([[-i,i],[0,0],[0,0]],int)
        else:
            ranges = nu.array([[-i+1,i],[0,0],[0,0]],int)
        return ranges
    
    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        rn=nu.zeros_like(r)
        x, y = r[0], r[1]
        rad = nu.sqrt(x**2+y**2)
        phi = phival(x,y)
        rn[0] = rad * cos( phi+n[0]*self.angle )
        rn[1] = rad * sin( phi+n[0]*self.angle )
        rn[2] = r[2]
        return rn

    def tensor(self,r,n):
        """ Dyadic tensor at r for symmetry operation n."""
        rn = self.transform(r,n)
        rad = nu.sqrt(r[0]**2+r[1]**2)
        assert rad>1E-10
        T = nu.array([[r[0]*rn[0]+r[1]*rn[1],-(r[0]*rn[1]-r[1]*rn[0]),0],\
                      [r[0]*rn[1]-r[1]*rn[0],  r[0]*rn[0]+r[1]*rn[1] ,0],\
                      [       0,                      0,         rad**2]])/rad**2
        return T   
    
    def rotation_of_axes(self,n):
        """ Rotation of axes for symmetry operation n."""
        angle = n[0]*self.angle
        R = nu.array([[cos(angle),sin(angle),0],[-sin(angle),cos(angle),0],[0,0,1]])
        return R  


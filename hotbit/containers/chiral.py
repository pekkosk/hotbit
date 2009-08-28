import numpy as nu
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
        assert type==self.type
        self.atoms = proxy(atoms)
        self.angle = None 
        self.height = atoms.get_cell()[2,2] 
        
    def __repr__(self):
        if self.angle==None:
            raise AssertionError('Chiral angle was not set yet.')
        if self.angle<1E-14:
            x='Chiral: angle=0.0, height=%.4f Ang' %self.height
        else:
            x='Chiral: angle=%.4f (2*pi/%.2f), height=%.4f Ang' %(self.angle,2*nu.pi/self.angle,self.height)
        return x

        
    def set(self,angle=None,height=None,scale_atoms=False,container=None):
        """
        Reset angle or height, and maybe scale atoms.
         
        @param: height   Height of the primitive cell in z-direction
        @param: angle    angle (in radians) of rotation
        """
        if container!=None:
            # copy container
            assert angle==None and height==None and scale_atoms==False
            self.set(angle=container.angle,height=container.height)
        else:
            if not scale_atoms:
                if angle!=None: self.angle = angle
                if height!=None: self.height = height
            else:
                if angle is None:
                    da = 0.0
                else:
                    if self.angle==None:
                        raise AssertionError('Positions cannot be scaled; initial angle was not given.')
                    da = angle - self.angle
                    self.angle = angle
                    
                old_height = self.height
                if height != None:
                    self.height = height
                    
                newr = []
                for r in self.atoms.get_positions():
                    x,y = r[0],r[1]
                    rad = nu.sqrt( x**2+y**2 )
                    frac = r[2]/old_height
                    # twist atoms z/h * da (more)
                    newphi = mix.phival(x,y) + frac * da
                    newr.append( [rad*nu.cos(newphi),rad*nu.sin(newphi),frac*self.height] )
                self.atoms.set_positions(newr)     
        
        self.atoms.set_pbc((False,False,True))
        self.atoms.set_cell(self.get_ase_cell())
        
    def __eq__(self,other):
        if isinstance(other,Chiral) and abs(self.angle-other.angle)<1E-12 \
           and abs(self.height-other.height)<1E-12:
            return True
        else:
            return False
    
    def get_ase_cell(self):
        """ cell used for visualization """
        if self.angle==None:
            raise AssertionError('Chiral angle is not set yet.')
        l = max(self.atoms.get_positions()[:,0])*1.5
        cell = nu.array( [[l,0,0],[l*cos(self.angle),l*sin(self.angle),0],[0,0,self.height]] )
        return cell
        
    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        if self.angle==None or self.height==None:
            raise RuntimeError("Chiral's angle or height is not set yet.")
        return nu.array( [[0,0],[0,0],[-nu.Inf,nu.Inf]] )
    
    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        rn=nu.zeros_like(r)
        x, y = r[0], r[1]
        rad = nu.sqrt(x**2+y**2)
        phi = phival(x,y)
        rn[0] = rad * cos( phi+n[2]*self.angle )
        rn[1] = rad * sin( phi+n[2]*self.angle )
        rn[2] = r[2] + n[2]*self.height
        return rn

    def tensor(self,r,n):
        """ Dyadic tensor at r for symmetry operation n."""
        rn = self.transform(r,n)
        rad = nu.sqrt(r[0]**2+r[1]**2)
        #assert rad>1E-10
        if rad>1E-10:
            T = nu.array([[r[0]*rn[0]+r[1]*rn[1],-(r[0]*rn[1]-r[1]*rn[0]),0],\
                          [r[0]*rn[1]-r[1]*rn[0],  r[0]*rn[0]+r[1]*rn[1] ,0],\
                          [       0,                      0,         rad**2]])/rad**2
        else:
            # atom on z-axis
            T = nu.zeros((3,3))
            T[2,2] = 1
        return T   
    
    def rotation_of_axes(self,n):
        """ Rotation of axes for symmetry operation n."""
        angle = n[2]*self.angle
        R = nu.array([[cos(angle),sin(angle),0],[-sin(angle),cos(angle),0],[0,0,1]])
        return R  
    

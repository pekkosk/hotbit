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
        self._set_table()
        
    def _set_table(self): 
        self.table = [{'M':1},{'M':1},{'M':nu.Inf}]
        
    def __repr__(self):
        if self.angle==None:
            raise AssertionError('Chiral angle was not set yet.')
        if self.angle<1E-14:
            x='Chiral: angle=0.0, height=%.4f Ang' %self.height
        else:
            x='Chiral: angle=%.4f (2*pi/%.2f), height=%.4f Ang' %(self.angle,2*nu.pi/self.angle,self.height)
        return x
    
    
    def get(self,key):
        """
        Return current angle or height
        """
        if key not in ['angle','height']:
            raise AssertionError('Invalid keyword %s' %key)
        if key=='angle':
            return self.angle
        elif key=='height':
            return self.height             
    

        
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
                    newphi = phival(x,y) + frac * da
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
        """ Rotate r by n2*angle + translate (in z) by n2*height."""
        R = self.rotation(n)
        return nu.dot(R,r) + (0,0,n[2]*self.height)
    
    def rotation(self,n):
        """ Return the (active) rotation matrix for symmetry operation n. """
        angle = n[2]*self.angle
        R = nu.array([[cos(angle),-sin(angle),0],[sin(angle),cos(angle),0],[0,0,1]])
        return R 
    

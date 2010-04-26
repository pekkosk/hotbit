import numpy as nu
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
        self.angle = None 
        self.height = atoms.get_cell()[2,2]
        self.x = None
        
    def __repr__(self):
        if self.angle==None:
            raise AssertionError('Chiral angle was not set yet.')
        if self.angle<1E-14:
            x='DoubleChiral: angle=0.0, height=%.4f Ang, x=%.4f' %(self.height,self.x)
        else:
            x='DoubleChiral: angle=%.4f (2*pi/%.2f), height=%.4f Ang, x=%.4f' %(self.angle,2*nu.pi/self.angle,self.height,self.x)
        return x
    
    
    def get(self,key):
        """
        Return current angle, height, or x
        """
        if key=='angle':
            return self.angle
        elif key=='height':
            return self.height             
        elif key=='x':
            return self.x
        else:
            raise AssertionError('Invalid keyword %s' %key)
    

        
    def set(self,angle=None,height=None,x=None,scale_atoms=False,container=None):
        """
        Reset angle, height, or x, and maybe scale atoms.
         
        parameters:
        ===========
        height:   Height of the primitive cell in z-direction
        angle:    angle (in radians) of rotation
        x:        fractional translation offset related to 180 rotation
        """
        if container!=None:
            # copy container
            assert angle==None and height==None and x==None and scale_atoms==False
            self.set(angle=container.angle,height=container.height,x=container.x)
        else:
            if not scale_atoms:
                if angle!=None: self.angle = angle
                if height!=None: self.height = height
                if x!=None: self.x = x
            else:
                if x!=None:
                    raise AssertionError('It is probably illegal to change x. This changes the symmetry, right?')
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
        
        self.atoms.set_pbc((True,False,True))
        self.atoms.set_cell(self.get_ase_cell())
        
    def __eq__(self,other):
        if isinstance(other,DoubleChiral) and abs(self.angle-other.angle)<1E-12 \
           and abs(self.height-other.height)<1E-12 \
           and abs(self.x-other.x)<1E-12:
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
        if self.angle==None or self.height==None or self.x==None:
            raise RuntimeError("Chiral's angle, height or x is not set yet.")
        return nu.array( [[0,1],[0,0],[-nu.Inf,nu.Inf]] )
    
    def transform(self,r,n):
        """ See init doc for transformation."""
        R = self.rotation(n)
        return nu.dot(R,r) + (0,0,(n[2]+n[0]*self.x)*self.height)
                         
                         
    def rotation(self,n):
        """ Return the (active) rotation matrix for symmetry operation n. """
        angle = n[2]*self.angle + n[0]*(nu.pi+self.angle*self.x)
        R = nu.array([[cos(angle),-sin(angle),0],[sin(angle),cos(angle),0],[0,0,1]])
        return R 
    

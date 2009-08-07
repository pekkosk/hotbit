from ase import Atoms as ase_Atoms
import numpy as nu
from box.mix import  phival
from math import sin,cos
from box import mix
from copy import copy 
from weakref import proxy


class Atoms(ase_Atoms):
    
    
    def __init__(self, symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 cell=None, pbc=None,
                 constraint=None,
                 calculator=None,
                 container='Bravais'):
        """ 
        Modified atoms class for hotbit.
        
        Input parameters as for ase.Atoms, except for keyword container.
        
        @param container: either container type ('Bravais','Wedge','Chiral'), or
                          dictionary describing the generalized unit cell, where                     
                          container['type'] selects the cell class; for other keywords,
                          look at the selected class
        """
        ase_Atoms.__init__(self,symbols=symbols,
             positions=positions, numbers=numbers,
             tags=tags, momenta=momenta, masses=masses,
             magmoms=magmoms, charges=charges,
             scaled_positions=scaled_positions,
             cell=cell, pbc=pbc,
             constraint=constraint,
             calculator=calculator)
        
        self._ranges = None                
        if type(container)==type(''):
            dict = {'type':container}
        else:
            dict = container.copy()
                
        # create the container instance
        assert 'type' in dict
        exec( 'self.container_class = %s' %dict['type'] )
        self.container = self.container_class(self,**dict)
        dict.pop('type')
        if dict!={}:
            self.container.set(**dict)
                    
        # these are just for shorter references
        self._transform = self.container.transform
        self._tensor = self.container.tensor
        self._rotation_of_axes = self.container.rotation_of_axes
                
        
    def set_container(self,**dict):
        '''
        Set the container class and its parameters
        
        @param dict: dictionary of container parameters
        '''
        if 'type' in dict:
            if dict['type']!=self.container.type: 
                raise AssertionError('Container type cannot be changed.')
        assert 'type' not in dict
        self.container.set(**dict)        
        
        
    def get_ranges(self):
        '''
        Return the ranges for symmetry operations in different directions.
        '''
        return self.container.get_ranges()
    
    
    def _check_symmetry_operation(self,n):
        '''
        Check that given symmetry operation is allowed.
        
        @param n: tuple for number of symmetry operations
        '''
        for i in range(3):
            a,b=self.container.ranges[i]
            if not a<=n[i]<=b:
                raise ValueError('Illegal symmetry operation: %i %i %i. For direction %i span [%i,%i] allowed.' %(n[0],n[1],n[2],i,a,b) )
            
            
    def transform(self,r,n):
        '''
        Transform position r according to symmetry operation n.
        
        @param r: position vector
        @param n: 3-tuple for symmetry operation.
        '''
        self._check_symmetry_operation(n)
        return self._transform(r,n)
        
        
    def tensor(self,r,n):
        '''
        Return the dyadic tensor 
        
                    d (R_j^n)_a
        T(jn)_ab = ------------- hat_a hat_b 
                    d (R_j)b
                    
        @param r: position vector
        @param n: symmetry operation 3-tuple
        '''
        self._check_symmetry_operation(n)
        return self._tensor(r,n)
    
    
    def rotation_of_axes(self,n):
        '''
        Return the 3x3 rotation matrix of coordination axes for given operation.
        
        @param n: 3-tuple for symmetry operation
        '''
        self._check_symmetry_operation(n)
        return self._rotation_of_axes(n)
                
                
    def extended_copy(self,n_list):
        """ Get copies of atoms for all listed symmetry operations n. 
        
        Return normal ase.Atoms -instance.
        """ 
        atoms2=None
        for n in n_list:
            self._check_symmetry_operation(n)
            atomsn=ase_Atoms(self)
            atomsn.set_positions( [self.transform(r,n) for r in self.get_positions()] )
            try:
                atoms2 += atomsn
            except:
                atoms2 = atomsn
        atoms2.set_pbc(False)
        atoms2.set_cell((1,1,1))
        return atoms2       


    def __eq__(self,other):       
        if ase_Atoms.__eq__(self,other):
            # for Bravais ase's Atoms.__eq__ is enough
            if self.container.type == 'Bravais':
                return True
            else:
                try:
                    same_container = self.same_container(other) 
                except:
                    raise ValueError('Comparing Bravais and non-Bravais containers should not happen. Check code.')
                if same_container:
                    return True
                else:
                    return False
        else:
            return False    
        
    def same_container(self,other):
        """ Check if atoms has the same container. """
        return self.container==other.container
        
    def copy(self):    
        """Return a copy."""
        cp = Atoms(container=self.container.type)
        cp += self
        cp.set_pbc( self.get_pbc() )
        cp.set_cell( self.get_cell() )
        if self.container.type!='Bravais':
            cp.set_container(container=self.container)
        return cp
        

class Bravais:
    
    def __init__(self,atoms,type):
        """
        Class for Bravais lattice, for hotbit.Atoms -class
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        """
        self.type = 'Bravais'
        assert type==self.type
        
        ranges = []
        for i in range(3):
            if atoms.get_pbc()[i]:
                ranges.append([-nu.Inf,nu.Inf])
            else:
                ranges.append([0,0])
        self.ranges = nu.array(ranges)    
        self.atoms = proxy(atoms)
        
    def set(self,**args):
        #print args
        #print args['container'].type
        #print args['container'].type!='Bravais'
        #print args
        #print 'container' in args and args['container'].type=='Bravais'
        if args!={}: # or not ('container' in args and args['container'].type=='Bravais'):
            raise NotImplementedError('For Bravais use set_pbc and set_cell normally.')
    
    def get_pbc(self):
        """ Return atoms's pbc as normal."""
        return self.atoms.get_pbc()
    
    def get_ase_cell(self):
        """ cell used for visualization """
        return self.atoms.get_cell()
    
    def __eq__(self,other):
        if isinstance(other,Bravais) and self.get_pbc()==other.get_pbc() \
           and self.get_cell()==other.get_cell():
            return True
        else:
            return False        
        
    def get_ranges(self):
        return self.ranges.copy()
    
    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        rn=r.copy()
        for a in range(3):
            rn = rn + n[a]*self.atoms._cell[a,:]
        return rn
        
    def tensor(self,r,n):
        """ Dyadic tensor at r for symmetry operation n."""
        return nu.eye(3)
            
    def rotation_of_axes(self,n):
        """ Rotation of axes for symmetry operation n."""
        return nu.eye(3)
    
    


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
        self.height = 1.0
        self.angle = 1E-12
        self.physical = True
                
            
            
    def set(self, angle=None, height=None, M=None, physical_angle=True, pbcz=None, scale_atoms=False, container=None):
        """ Only height can be reset, not angle. 

        @param: angle    angle (in radians) of the wedge (and M=None)
        @param: height   Height of the primitive cell in z-direction
        @param: M        set angle to 2*pi/M (with angle=None)
        @param: physical_angle (only for M=None) if angle is small, it does not be
                         exactly 2*pi/integer, i.e. situation has no physical meaning
                         (use for calculating stuff continuously)
        @param: pbcz     True if wedge is periodic in z-direction
        @param: scale_atoms Scale atoms according to changes in parameters 
        """
        if container!=None:
            assert angle==None and height==None and M==None and pbcz==None
            self.set(angle=container.angle,height=container.height,\
                     physical_angle=container.physical, pbcz=container.atoms.get_pbc()[2])
        
        if angle!=None or M!=None:
            assert not scale_atoms
            assert not (angle!=None and M!=None)
            if angle is None:
                assert isinstance(M,int)
                self.angle = 2*nu.pi/M
                self.M = M
            else:
                self.M = int( round(2*nu.pi/angle) )
                self.angle = angle
                
            # check parameters
            self.physical = physical_angle
            if self.angle<1E-6:
                raise Warning('Too small angle (%f) may bring rounding problems.' %self.angle)
            if self.angle>nu.pi:
                raise AssertionError('angle>pi')
            if nu.abs(self.M-2*nu.pi/self.angle)>1E-12 and physical_angle: 
                raise AssertionError('angle not physical; angle != 2*pi/M')
            if not physical_angle and self.M<20:
                raise AssertionError('Too large, non-physical angle.')
    
            i = self.M/2
            if nu.mod(self.M,2)==1:
                self.ranges = nu.array([[-i,i],[0,0],[0,0]],int)
            else:
                self.ranges = nu.array([[-i+1,i],[0,0],[0,0]],int)
                
        if height!=None:
            if scale_atoms:
                r = self.atoms.get_positions()
                r[:,2] = r[:,2] * height/self.height
                atoms.set_positions(r)
            self.height = height
            self.atoms.set_cell( self.get_ase_cell() )
            
        if pbcz!=None:
            self.atoms.set_pbc((True,False,pbcz))
            
        self.atoms.set_cell(self.get_ase_cell())  
                   
            
    def __eq__(self,other):
        if isinstance(other,Wedge) and abs(self.height-other.height)<1E-12 \
           and abs(self.angle-other.angle)<1E-12 and self.atoms.get_pbc()==other.atoms.get_pbc():
            return True
        else:
            return False
    
    def get_ase_cell(self):
        """ cell used for visualization """
        l = max(self.atoms.get_positions()[:,0])*1.5
        return nu.array( [[l,0,0],[l*cos(self.angle),l*sin(self.angle),0],[0,0,self.height]] )
        
    def get_ranges(self):
        return self.ranges.copy()
    
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
        self.angle = 0.0
        self.height = 1.0
        self.ranges = nu.array( [[0,0],[0,0],[-nu.Inf,nu.Inf]] )
        
        
    def set(self,angle=None,height=None,scale_atoms=False,container=None):
        """
        Reset angle or height, and maybe scale atoms.
         
        @param: height   Height of the primitive cell in z-direction
        @param: angle    angle (in radians) of rotation
        """
        if container!=None:
            self.set(angle=container.angle,height=container.height)
        if not scale_atoms:
            if angle!=None: self.angle = angle
            if height!=None: self.height = height
        else:
            if angle is None:
                da = 0.0
            else:
                da = angle - self.angle
                self.angle = angle
                
            if height is None:
                fact = 1.0
            else:
                fact = height/self.height
                self.height = height
                
            newr = []
            for r in atoms.get_positions():
                x,y = r[0],r[1]
                rad = nu.sqrt( x**2+y**2 )
                r[2] = fact*r[2]
                # twist atoms z/h * da *more*
                newphi = mix.phival(x,y) + r[2]/self.height * da
                newr.append( [rad*nu.cos(newphi),rad*nu.sin(newphi),r[2]] )
            atoms.set_positions(newr)     
        
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
        l = max(self.atoms.get_positions()[:,0])*1.5
        return nu.array( [[l,0,0],[l*cos(self.angle),l*sin(self.angle),0],[0,0,self.height]] )
        
    def get_ranges(self):
        return self.ranges.copy()
    
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
        assert rad>1E-10
        T = nu.array([[r[0]*rn[0]+r[1]*rn[1],-(r[0]*rn[1]-r[1]*rn[0]),0],\
                      [r[0]*rn[1]-r[1]*rn[0],  r[0]*rn[0]+r[1]*rn[1] ,0],\
                      [       0,                      0,         rad**2]])/rad**2
        return T   
    
    def rotation_of_axes(self,n):
        """ Rotation of axes for symmetry operation n."""
        angle = n[2]*self.angle
        R = nu.array([[cos(angle),sin(angle),0],[-sin(angle),cos(angle),0],[0,0,1]])
        return R  
    
                  
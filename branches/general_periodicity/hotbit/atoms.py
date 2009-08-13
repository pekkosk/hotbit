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
        r = self.container.get_ranges()
        for i in range(3):
            a,b=r[i]
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
                
                
    def extended_copy(self,n):
        """ Get copies of atoms for all listed symmetry operations n.
        
        @param: n  can be a list of 3-tuples for transformations, or 3-tuple
                   for telling how many copies in each direction is made.  
        
        Return normal ase.Atoms -instance.
        """ 
        r = self.container.get_ranges()
        if isinstance(n,list):
            n_list = n.copy()
        elif isinstance(n,tuple):
            a = []
            for i in range(3):
                if r[i,0]==-nu.inf:
                    a.append(0)
                else:
                    M = r[i,1] + 1
                    # try to start copies from primitive cell 0 first
                    if n[i]>M:
                        M = r[i,1] - r[i,0] + 1
                        a.append(r[i,0])
                    else:
                        a.append(0)
                    #M = r[i,1] - r[i,0] + 1
                    if n[i]>M:
                        raise AssertionError('Too many extended copies for direction %i.' %i)
                            
            n_list=[]
            for n1 in range(n[0]):
                for n2 in range(n[1]):
                    for n3 in range(n[2]):
                        n_list.append( (a[0]+n1,a[1]+n2,a[2]+n3) )
        
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
        #print self._cell
        #print other._cell
        #print ase_Atoms.__eq__(self,other)
        if ase_Atoms.__eq__(self,other):
            # for Bravais ase's Atoms.__eq__ is enough
            if self.container.type == 'Bravais':
                return True
            else:
                if hasattr(other,'container'):
                    return self.same_container(other)
                else:
                    raise AssertionError('Comparing Bravais and non-Bravais containers should not happen. Check the code.')
        else:
            return False    
        
    def same_container(self,other):
        """ Check if atoms has the same container. """
        return self.container==other.container
        
    def copy(self):    
        """Return a copy."""
        cp = Atoms(container=self.container.type)
        cp += self
        # set cell and pbc for initialization
        cp.set_pbc( self.get_pbc() )
        cp.set_cell( self.get_cell() )
        if self.container.type!='Bravais':
            cp.set_container(container=self.container)
        # reset cell (for ase and visualization use) exactly the same
        # (ase cell in set_container is taken from present atom positions,
        # even though originally it might have been set earlier)
        assert all( self.get_pbc()==cp.get_pbc() ) 
        cp.set_cell( self.get_cell() )
        return cp
        

class Bravais:
    
    def __init__(self,atoms,type):
        """
        Class for Bravais lattice, for hotbit.Atoms -class
        
        More documentation for the methods can be found from hotbit.Atoms -class. 
        """
        self.type = 'Bravais'
        assert type==self.type
        self.atoms = proxy(atoms)
        
    def __repr__(self):
        pbc=self.atoms.get_pbc()
        cell=self.atoms.get_cell()
        d = []
        for a in range(3):
            d.append( nu.linalg.norm(cell[a,:]) )
            
        a12 = nu.dot(cell[0],cell[1])/(d[0]*d[1])
        a13 = nu.dot(cell[0],cell[2])/(d[0]*d[2])
        a23 = nu.dot(cell[1],cell[2])/(d[1]*d[2])
        x='Bravais: pbc:[%i,%i,%i], ' %(pbc[0],pbc[1],pbc[2])
        x+='cell:[%.2f,%.2f,%.2f] Ang, angles(12,13,23):[%.2f,%.2f,%.2f]' %(d[0],d[1],d[2],a12,a13,a23) 
        return x
        
    def set(self,**args):
        if args!={}: 
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
        """ Return ranges for symmetry operations. """
        ranges = []
        for i in range(3):
            if self.atoms.get_pbc()[i]:
                ranges.append([-nu.Inf,nu.Inf])
            else:
                ranges.append([0,0])
        return nu.array(ranges)
    
    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        rn=r.copy()
        cell = self.atoms.get_cell()
        for a in range(3):
            rn = rn + n[a]*cell[a,:]
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
        self.height = None
        self.angle = None
        self.physical = None
        
    def __repr__(self):
        x='Wedge: angle=%.4f (2*pi/%.2f, ' %(self.angle,2*nu.pi/self.angle)
        if self.physical:
            x+='physical), '
        else:
            x+='not physical), '
        x+='height=%.4f Ang ' %self.height
        if self.atoms.get_pbc()[2]:
            x+='(pbc)'
        else:
            x+='(no pbc)'
        return x
            
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
            self.physical = physical_angle
            if self.angle<1E-6:
                raise Warning('Too small angle (%f) may bring rounding problems.' %self.angle)
            if self.angle>nu.pi:
                raise AssertionError('angle>pi')
            if nu.abs(self.M-2*nu.pi/self.angle)>1E-12 and physical_angle: 
                raise AssertionError('angle not physical; angle != 2*pi/M')
            if not physical_angle and self.M<20:
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
        cell = atoms.get_cell()
        self.height = cell[2,2] 
        
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
        
    def get_ranges(self):
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
    
                  
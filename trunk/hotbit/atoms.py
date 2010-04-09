from ase import Atoms as ase_Atoms
import numpy as nu
from box import mix
from copy import copy 
from hotbit.containers import *


class Atoms(ase_Atoms):
    
    
    def __init__(self, symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 cell=None, pbc=None,
                 constraint=None,
                 calculator=None,
                 container='Bravais',
                 atoms=None):
        """ 
        Modified atoms class for hotbit.
        
        Input parameters as for ase.Atoms, except for keyword container.
        
        @param container: either container type ('Bravais','Wedge','Chiral' or any other
                          container from hotbit/containers/*.py), or
                          dictionary describing the generalized unit cell, where                     
                          container['type'] selects the cell class; for other keywords,
                          look at the selected classes
        """
        ase_Atoms.__init__(self,symbols=symbols,
             positions=positions, numbers=numbers,
             tags=tags, momenta=momenta, masses=masses,
             magmoms=magmoms, charges=charges,
             scaled_positions=scaled_positions,
             cell=cell, pbc=pbc,
             constraint=constraint,
             calculator=calculator)
        
        if type(atoms)!=type(None):
            self += atoms
            self.set_pbc( atoms.get_pbc() )
            self.set_cell( atoms.get_cell() )
                        
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
        self._rotation = self.container.rotation
                
        
    def set_container(self,**cont):
        '''
        Set the container class and its parameters
        
        @param dict: dictionary of container parameters
        '''
        if 'type' in cont:
            if cont['type']!=self.container.type: 
                raise AssertionError('Container type cannot be changed.')
        assert 'type' not in cont
        self.container.set(**cont)        
        
        
    def get_symmetry_operation_ranges(self):
        '''
        Return the ranges for symmetry operations in different directions.
        '''
        return self.container.get_symmetry_operation_ranges()
    
    
    def _check_symmetry_operation(self,n):
        '''
        Check that given symmetry operation is allowed.
        
        @param n: tuple for number of symmetry operations
        '''
        r = self.container.get_symmetry_operation_ranges()
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
        assert False
        self._check_symmetry_operation(n)
        return self._tensor(r,n)
    
    
    def rotation_of_axes(self,n):
        '''
        Return the 3x3 rotation matrix of coordination axes for given operation.
        
        @param n: 3-tuple for symmetry operation
        '''
        assert False
        self._check_symmetry_operation(n)
        return self._rotation_of_axes(n)
    
    def rotation(self,n):
        self._check_symmetry_operation(n)
        return self._rotation(n)
                
                
    def extended_copy(self,n):
        """ Get copies of atoms for all listed symmetry operations n.
        
        @param: n  i) list of 3-tuples for transformations explicitly
                      e.g. [(0,0,0),(1,0,0),(0,1,0)]
                   ii) 3-tuple for the number of transformations 
                       (for infinite extensions -> start from zero;
                        for finite extensions -> symmetry ops around zero)
                       e.g. (10,10,1)
                   iii) 3-tuple of 2-tuples that give the ranges for symmetry ops 
                       e.g. ((-3,3),(-3,3),1)
        
        Return normal ase.Atoms -instance.
        """ 
        r = self.container.get_symmetry_operation_ranges()
        if isinstance(n,list):
            n_list = copy(n)
        elif isinstance(n,tuple):
            ops = []
            for i in range(3):
                if isinstance(n[i],tuple):
                    ops.append( nu.arange(n[i][0],n[i][1]+1) )
                elif isinstance(n[i],int):
                    assert n[i]>0
                    rng = nu.arange(0,n[i])
                    if r[i,0]!=-nu.Inf:
                        # if nr. of ops is finite, center copies wrt. 0
                        N = r[i,1] - r[i,0] + 1
                        if len(rng)>N:
                            raise AssertionError('Too many extended copies (%i) for direction %i. Only %i is allowed.' %(len(rng),i,N) )
                        rng = rng - len(rng)/2
                    ops.append( rng )
                            
            n_list=[]
            for n1 in ops[0]:
                for n2 in ops[1]:
                    for n3 in ops[2]:
                        n_list.append( (n1,n2,n3) )
        
        atoms2=None
        for n in n_list:
            self._check_symmetry_operation(n)
            atomsn = ase_Atoms()
            atomsn += self
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
        assert nu.all( self.get_pbc()==cp.get_pbc() ) 
        cp.set_cell( self.get_cell() )
        return cp
        
    def __imul__(self, m):
        raise NotImplementedError('*= not implemented yet, use extended_copy.')



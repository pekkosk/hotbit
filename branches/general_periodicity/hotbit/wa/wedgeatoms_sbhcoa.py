#
# safer-but-has-copies-of-ase pices!!!
#
# Module under development!

from math import pi as _pi
 
import numpy as npy

from ase import Atoms
from ase.data import atomic_numbers, chemical_symbols, atomic_masses

class WedgeAtoms(Atoms):
    # TODO: not tested!
    #
    # TODO: rewrite all system Atoms' methods
    # TODO: check what other  Atoms' methods needs no be rewritten?
    #
    # DOC: add ASCII picture of wedge geometry
    #
    # TODO: angle in degrees?
    # TODO: angle -> wedge_angle
    # FIXME: to remove pbc attribute is a bad idea
    # find a good solution for this, pbc in x and y 
    #
    # 
    # TODO: shell angle be a num array?! deepcopy?!
    #
    # 
    # TODO: implement __slots__ 
    # __slots__ makes Atoms consume less memory!
        
    """
    Wedge Atoms Object.
        
    !! safer-but-has-copies-of-ase pices!! 
    Ggain: has not pbc attribute, no possibility to change it. 
    
    (But, if user of Wedge Atoms assumes that the object has
    an accessible pbc attribute, can be a source of problem.)
    NO BUT! -> proper usage of Wedge Atoms in all cases! 
    
     
     
        
    self.is_wedge_atoms_instance() == True
    self.is_wedge_atoms_instance == True ?!
    pbc attribute is [None, None, pbc_z]
    
    Attributes Introduced and with Shifted meaning: 
    
    angle:
    
    pbc: inherited from Atoms.pbc, BUT 1. 2. 
        
        1.Should not be manipulated direcly, 
        use set_pbc_z() get_pbc_z() instead. See also Comment 1. below  
        
        2. Here it is different! The Periodic Boundary Conditions 
        are represented by ints (provided that _pbc_repr_int_not_bool 
        == True, it is set so by __init__()): 
        
        *New representation:
        
        [2,2, pbc_z],
        
        where 'pbc_z' is 1 or 0, meaning True or False for 
                  periodic boundary conditions along Z axis 
              '2' reflects the fast that  periodic boundary conditions 
                  along X and Y axis are always True by definition
                  of clas.
              
        e.g. [2,2,1] or  [2,2,0],
        See also __init__(), set_pbc_z(), get_pbc_z(),
        _pbc_repr_int_not_bool.    
    
    
    Example: 
    ...
    
    Comments:
    
    1. Example of what you should NOT do with pbc attribute. 
    As in ASE,
    
    >>> a = Atoms(molecule('CO'))
    >>> a.pbc = True  # Ok,  but
    >>> print a.get_pbc() # Fails
    
    the pbc attribute is not supposed to be changed directly:
     
    >>> w = WedgeAtoms(_pi / 6, molecule('CO'), pbc=True)
    >>> w.pbc = False # # Ok,  but
    >>> print w.get_pbc() # results in error 
    
    The best solution would be to add in ase.Atoms pbc as property!
    The lines 
    (self.pbc == other.pbc).all())  from __eq__(self, other):
    s += 'pbc=%s, ' % self.pbc.tolist()  from ase.Atoms.__repr__()
    Should be changed. 
    At the 3.0.0. implementation the are sources of problems.
    Can we rewrite it?
     
    2.__repr__() behave as Atoms.__repr__():  
    # Atoms.__repr__()
    # Does not work with b = eval(repr(a))!
    # http://docs.python.org/ says: this function [repr(object)] makes an attempt to return a string that would yield an object with the same value when passed to eval(), otherwise the representation is a string enclosed in angle brackets that contains the name of the type of the object together with additional information often including the name and address of the object. 
    >>> a = Atoms(molecule('CO'), pbc=True)
    >>> b = eval(repr(a)) # results in SyntaxError due to 'positions=...'
    >>> print v == w
    
     
    
    3.<no more comments>
    """
    
    # Static vars
    _min_angle = 0.
    # FIXME: what _max_angle should be? 
    _max_angle = _pi / 2
    
    def __init__(self,
                 # New argument
                 angle,
                 symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 # TODO: cell might need to be changed
                 cell=None,
                 pbc_z=None,
                 constraint=None,
                 calculator=None):
        
        """ 
        Sets _pbc_repr_int_not_bool to True.
        
        Parameters:
        
        @param angle: float, must be in range (_min_angle, _max_angle]
            wedge angle in radians
        @param symbols: same as in ase.Atoms.__init__()
        @param positions: same as in ase.Atoms.__init__()
        @param numbers: same as in ase.Atoms.__init__()
        @param tags: same as in ase.Atoms.__init__()
        @param momenta: same as in ase.Atoms.__init__()
        @param masses: same as in ase.Atoms.__init__()
        @param magmoms: same as in ase.Atoms.__init__()
        @param charges: same as in ase.Atoms.__init__()
        @param scaled_positions: same as in ase.Atoms.__init__()
        @param cell: same as in ase.Atoms.__init__()
        @param pbc_z: (exactly one)  bool (Default: False)
            Periodic Boundary Conditions flag along Z axis.
            Examples: True, False, 0, 1. 
        @param constraint: same as in ase.Atoms.__init__()
        @param calculator: same as in ase.Atoms.__init__()
        
        See  also class' doc string.
        """
        
        # TODO: why this is not working?
        # TODO: check if this work correct
        # should hide self.__dic__
        #__slots__ = ['pbc']
        #__slots__ = []
        
        self._is_wedge_atoms_instance = True
        self.set_angle(angle)
        #self.angle = angle # use is angle is a property 
        
        # DOC: add help
        self._pbc_repr_int_not_bool = True    
        
        if pbc_z is None:
            pbc_z = False 
        
        
        Atoms.__init__(self, symbols=symbols,
                 positions=positions, numbers=numbers,
                 tags=tags, momenta=momenta, masses=masses,
                 magmoms=magmoms, charges=charges,
                 scaled_positions=scaled_positions,
                 # TODO: cell might need to be changed
                 cell=cell,
                 pbc=pbc_z,
                 constraint=constraint,
                 calculator=calculator)
     
    def set_pbc(self, pbc):
        """
        Use set_pbc_z() instead to set pbc along Z axis.
        
        For backward compatibility with ase.Atoms and 
        should not be used unless you are not satisfied with set_pbc_z()
        
        (Example of not a proper usage:
        >>> w = WedgeAtoms(0.1, molecule('CO'), pbc=True)
        >>> w.set_pbc([False,False,False]) # is an error in)
        """
        
        # None will harm here if _pbc_repr_int_not_bool is false
        pbc_x_y_always_true = 2   

        if isinstance(pbc, (list, tuple)):
            # FIXME: exception
            #print 'in set_pbc'
            #print 'pbc:', pbc
            #print 'pbc.type:', type(pbc)
            #print 'pbc.tolist:', pbc.tolist()
            raise Exception
        
        if isinstance(pbc, int):
            pbc = (pbc_x_y_always_true, pbc_x_y_always_true, pbc)
        
        # otherwise it must be an instance of isinstance(pbc, (npy.ndarray,)):
        # which is ok  
        # or a bad user input 
        
#===============================================================================
#        if isinstance(pbc, (npy.ndarray,))
#            # FIXME: hardcoded data!
#            pbc = True
#            pbc = (pbc_x_y_always_true, pbc_x_y_always_true, pbc)
#===============================================================================

        if self._pbc_repr_int_not_bool:
            self._pbc = npy.array(pbc, int)
        else:
            self._pbc = npy.array(pbc, bool)
       
    def get_pbc(self):
        """
        Use get_pbc_z() instead to set pbc along Z axis.
        
        For backward compatibility with ase.Atoms and 
        should not be used unless you are not satisfied with get_pbc_z()
        
        Get periodic boundary condition flags (in internal format)
        """
        # FIXME: generate warning instead of print
        print 'Warning: %s.get_pbc() was used!' % self.__class__.__name__
        
        return self._pbc.copy()
        #return Atoms.get_pbc(self)
    
    def set_pbc_z(self, pbc_z):
        """
        Set Periodic Boundary Conditions flag along Z axis.
        @param pbc_z: (exactly one) bool
        """
        # The action is redirected to set_pbc(self, pbc)
        self.set_pbc(pbc_z)
    
    def get_pbc_z(self):
        """
        Get periodic boundary condition flag in z direction! 
        (in internal format)
        """
        
        #print self._pbc[0]
        #print self._pbc[1]
        #print bool(self._pbc[2])
        return bool(self._pbc[2]) 
    
    # DOC: this
    pbc_z = property(fget=get_pbc_z, fset=set_pbc_z, doc='periodic  boundary condition flag in z direction. Uses get/set methods')
    #pbc = property(fget=get_pbc_z, fset=set_pbc_z, doc='New pbc attribute!')
     
    def copy(self):
        """
        Return a copy.
        
        Note, implementation was partially copied from ase 3.0.0,
        from ase.Atoms.copy(self) in particular, 
        and then this-class-specific functionality was added.
        Thus, should be updated if ase changes! 
        """
        # TODOed: not tested -> seam to work 
        
        #wedgeatoms = WedgeAtoms(angle=self.get_angle(),
        #    # FIXME: pbc=self.get_pbc_z() ?! 
        #                        pbc_z=self.get_pbc_z(),
        #   # deep copy of all arrays, constraints, adsorbate_info...
        #    # must be done by ?!
        #    # FIXME: are atoms here ok??
        #                        symbols=Atoms.copy(self))
        
        # Copied from ase.Atoms
        import copy
        wedgeatoms = WedgeAtoms(angle=self.get_angle(), 
                                pbc_z=self.get_pbc_z(),
                                cell=self.cell)

        wedgeatoms.arrays = {}
        for name, a in self.arrays.items():
            wedgeatoms.arrays[name] = a.copy()
        wedgeatoms.constraints = copy.deepcopy(self.constraints)
        wedgeatoms.adsorbate_info = copy.deepcopy(self.adsorbate_info)
        return wedgeatoms
    
    def __repr__(self):      
        """
        Return informal object representation 
        (see Comments (1.) in the class docstring )
        
        Note, implementation was partially copied from ase 3.0.0,
        from ase.Atoms.__repr__(self) in particular, 
        and then this-class-specific functionality was added.
        Thus, should be updated if ase changes! 
        """
        
        # Copied from ase.Atoms.__repr__(self)
        num = self.get_atomic_numbers()
        N = len(num)
        if N == 0:
            symbols = ''
        elif N <= 60:
            # Distinct atomic numbers in num:
            dis = npy.concatenate(([0], npy.arange(1, N)[num[1:] != num[:-1]]))
            repeat = npy.append(dis[1:], N) - dis
            symbols = ''.join([chemical_symbols[num[d]] + str(r) * (r != 1)
                               for r, d in zip(repeat, dis)])
        else:
            symbols = ''.join([chemical_symbols[Z] for Z in num[:15]]) + '...'
        s = "Atoms(symbols='%s', " % symbols
        for name in self.arrays:
            if name == 'numbers':
                continue
            s += '%s=..., ' % name
        if abs(self.cell - npy.diag(self.cell.diagonal())).sum() < 1e-12:
            s += 'cell=%s, ' % self.cell.diagonal().tolist()
        else:
            s += 'cell=%s, ' % self.cell.tolist()            
        #s += 'pbc=%s, ' % self.pbc.tolist()
        if len(self.constraints) == 1:
            s += 'constraint=%s, ' % repr(self.constraints[0])
        if len(self.constraints) > 1:
            s += 'constraint=%s, ' % repr(self.constraints)
        if self.calc is not None:
            s += 'calculator=%s(...), ' % self.calc.__class__.__name__
        
        s = s[:-2] + ')'
    
        # new
        # this behavior relies that s 
        # 1. starts from 'Atoms(' and 
        # 2. ends with ')'     
        
        # add class name in front
        # insert  'angle = ...'  after 'Atoms(' at the begining  
        s = self.__class__.__name__+ '(angle = %s, ' % repr(self.get_angle()) + s[6:]
        # insert 'pbc_z= ...'
        s = s[:-1] + 'pbc_z='+ repr(self.get_pbc_z()) + ')'
        return s
    
#===============================================================================
#    def __repr__(self):
#        """
#        Return informal object representation 
#        (see Comments (1.) in the class docstring )
#        
#        Note, accurate behavior relies ase.Atoms.__repr__(self) return, in particular,
#        1. it starts from 'Atoms(' and 2. it contains pattern 'pbc=[ <zero or more of any character >]'     
#        """
#        # FIXME: accurate behavior relies ase.Atoms.__repr__(self) return  
#        
#        s = Atoms.__repr__(self)
#        
#        import re
#        pattern = r'pbc=\[.*\]' # pbc=[ <zero or more of any character >]
#        if not (s.startswith('Atoms(') and re.search(pattern, s)):
#            # FIXME: generate warning 
#            s += '*Warning: this might be not accurate! If ase.Atoms.__repr__() changes,'+ self.__class__.__name__+'.__repr__() might need changes' 
#        # add class name in front
#        # insert  'angle = ...'  after 'Atoms(' at the begining  
#        s = self.__class__.__name__+ '(angle = %s, ' % repr(self.get_angle()) + s[6:]
#        # insert 'pbc_z= ...'
#        s = re.sub(pattern,'pbc_z='+ repr(self.get_pbc_z()), s)
#        return s
#===============================================================================
    
    def __eq__(self, other):
        """
        Check for identity of two wedge atoms  objects.

        Identity means: same wedge angles, positions, atomic numbers, unit cell and
        periodic boundary conditions.
        
        Note, implementation was partially copied from ase 3.0.0,
        from ase.Atoms.__eq__(self, other) in particular, 
        and then this-class-specific functionality was added.
        Thus, should be updated if ase changes! 
        """        
        
        # Mosly copied from ase.Atoms
        a = self.arrays
        b = other.arrays
        return (len(self) == len(other) and
                (a['positions'] == b['positions']).all() and
                (a['numbers'] == b['numbers']).all() and
                (self.cell == other.cell).all() and
                #this line chanhed 
                (self.get_pbc_z() == other.get_pbc_z()) and
                #this line added only
                self.get_angle() == other.get_angle())
        
        #return (Atoms.__eq__(self, other) and 
        #        self.get_angle() == other.get_angle()) 
        
        
    
    def get_angle(self):
        """
        Return Wedge angle in radians
        """
        return self._angle
    
    def set_angle(self, angle):
        """
        Set Wedge angle 
        @param angle: in radians
        """
        # FIXME: add type check!
        # FIXME: should it be a set-able ?!
        # print 'in set_angle()'
        if self._min_angle < angle and angle <= self._max_angle:  
            if self.are_atoms_in_wedge():
                self._angle = angle
            else: 
                # FIXME: Exception
                raise Exception2
        else:
            # FIXME: Exception
            print 'angle:', angle
            raise Exception1
        
    # Should it be a property ?!
    # with __slot__ there is not much point:
    #angle = property(fget=get_angle, fset=set_angle, doc='Wedge angle in radians')
        
    
    def is_wedge_atoms_instance(self):
        """
        Return True if self is WedgeAtoms instance 
        """
        return self._is_wedge_atoms_instance
        #return self._wedge_atoms_flag
    
    # isn't function not enough?
    # FIXME: WTF? _is_wedge_atoms_instance not a hidden argument!
    is_wedge_atoms_instance = property(fget=is_wedge_atoms_instance, doc='True if self is WedgeAtoms instance')
    
    
#===============================================================================
#    def write_angle(self):
#        # FIXME: use other output, __repr__ of angle
#        """
#        Writes Wedge angle to stdout
#        """
#        if self.is_wedge_atoms_instance: 
#            print 'Wedge angle = ', self.angle
#        else: 
#            print 'This is not a Wedge Atoms instance'   
#===============================================================================
    
    
    def are_atoms_in_wedge(self):
        """ 
        Return bol
        Checks weather all atoms are located
        within the wedge
        """
        # FIXME: implement
        # FIXME: add whenever the coordinates are update
        # FIXME: add property for this function value?
        #print 'in are_atoms_in_wedge()'
        return True

    
    def get_vector(self, i, j):
        # FIXME: implement
        """
        Returns minimun image vector i,j
        See: elements.vector
        """
        return None

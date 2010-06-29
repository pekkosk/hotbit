# Module under development!

from numpy import array as numpy_array
from numpy import zeros as numpy_zeros 
from numpy.linalg import norm as numpy_linalg_norm
from numpy import nan as numpy_nan
from numpy import array as numpy_array  
from numpy import ndarray as numpy_ndarray
from numpy import Inf as numpy_Inf

from numpy import float64 as numpy_float64 

from ase import Atoms
#from ase.data import atomic_numbers, chemical_symbols, atomic_masses

from math import pi, sqrt
from math import sin as sin
from math import cos as cos

class WedgeAtoms(Atoms):
    # TODO: not all fetures tested!
    # DOC: check spelling of doc strings
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
    # TODO: add the from_angle  
    # from_angle=0 # New argument
    #    @param from_angle: Wedge is located at  
    #        phi belongs to (from_angle, from_angle + angle),
    #        where phi is a polar angle of
    #       polar coordinates (x,y) --> (r, phi)
        
        
    """
    Wedge Atoms Object.
        
    self.is_wedge_atoms_instance == True !
    pbc attribute is [None, None, pbc_z]
    
    Some of inherited methods may not have good sense
    Some of inherited methods may/will/should result in error
    
    Attributes Introduced and with Shifted meaning: 
    
    angle:
    
    pbc: inherited from Atoms.pbc, BUT 1. 2. 
        
        1.Should not be manipulated directly, 
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
    
    cell: inherited from Atoms.cell, can be initialized 
        and manipulated as in Atoms, BUT it does NOT carry much 
        meanings in Wedge geometry.
    
    
    Example: 
    ...
    
    Comments:
    
    1.and 2. - pbc,cell attributes are not properties!
     
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
     
    2. Same applies to cell attribute:
    
    >>> w = WedgeAtoms(_pi / 6, molecule('CO'), pbc=True)
    >>> w.cell=[2,2,2]
    >>> print w.get_cell() # results in error
    
    3.__repr__() behave as Atoms.__repr__():  
    # Atoms.__repr__()
    # Does not work with b = eval(repr(a))!
    # http://docs.python.org/ says: this function [repr(object)] makes 
    an attempt to return a string that would yield an object with the same 
    value when passed to eval(), otherwise the representation is a string 
    enclosed in angle brackets that contains the name of the type of the 
    object together with additional information often including the name 
    and address of the object. 
    
    >>> a = Atoms(molecule('CO'), pbc=True)
    >>> b = eval(repr(a)) # results in SyntaxError due to 'positions=...'
    >>> print v == w
    
     
    
    4.<no more comments>
    """
    
    # Static vars
    _min_angle = 0.
    # FIXME: what _max_angle better to use? 
    #_max_angle = pi / 2
    _max_angle = 2 * pi
    
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
                 # New argument
                 #pbc_z=False, 
                 # ->> Errors: a = Atoms(pbc = True); w = WedgeAtoms(radians(30), a) !!
                 pbc_z=None,
                 constraint=None,
                 calculator=None,
                 # New argument
                 number_of_cells=None
                 ):
        
        # DOC: from_angle
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
        self.angle = angle 
        
        # DOC: cell 
        #Format
        #cell = r_max
        #cell = [r_max]
        #cell = [r_max, delta_z]
        #cell = [r_max, [0,0,delta_z]]
        
        
        # FIXME: delta_z and r_max to be set together only
        
        # default
        delta_z = 1
        r_max = 1
        
        if (cell is not None) and (isinstance(cell, int)):
            cell = [cell]
         
        if (cell is not None) and (isinstance(cell, (list, tuple))):
            r_max = cell[0]
            if len(cell) > 1:
                if (isinstance(cell[1], int)):
                    delta_z = cell[1]
                elif (isinstance(cell[1], (list, tuple))):
                    if (cell[1][0] == 0) and (cell[1][1] == 0): #and (cell[1][2] != 0):
                        delta_z = cell[1][2]
                    else:
                        raise NotImplementedError('Such z-geometry '
                            + 'is not implemented yet')
                else:
                    raise TypeError('cell attribute should be'
                            + ' int or tuple/list')
                    
            cell = [[r_max, 0, 0],
                [r_max * cos(angle), r_max * sin(angle), 0],
                [0, 0, delta_z]]
        
        elif (cell is None):    
            cell = [[r_max, 0, 0],
                [r_max * cos(angle), r_max * sin(angle), 0],
                [0, 0, delta_z]]
        
        # FIXME: add checks to preserve direct manipulations
        # numpy_ndarray can be when in copying...
        elif not (isinstance(cell, numpy_ndarray)):
            #print type(cell) 
            raise TypeError('cell attribute should be'
                            + ' int or tuple/list')
        
        
        # FIXME: cell is output twice by __rep__!!
        
        
        # DOC: add help
        self._pbc_repr_int_not_bool = False
        #self._pbc_repr_int_not_bool = True    
        
        #if pbc_z is None:
        #    pbc_z = False 
        
        #FIXEME:
        # if isinstance(symbols, Atoms):
        #    atoms = symbols
        #    symbols = None    
        
        # If WedgeAtoms are not constructed from Atoms
        if  (pbc_z is None) and (not isinstance(symbols, Atoms)):
            pbc_z = False
        #    atoms = symbols
        #    symbols = None    
        
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
        
        # should go after set_pbc!
        self.set_number_of_cells(number_of_cells)
     
    def set_pbc(self, pbc):
        #DOC: update this!
        """
        The value of atoms.pbc protected!
        
        The pbc value of numpy type used incorrelty can led to 
        error!
        
        Use set_pbc_z() instead to set pbc along Z axis.
        
        For backward compatibility with ase.Atoms and 
        should not be used unless you are not satisfied with set_pbc_z()
        
        (Example of not a proper usage:
        >>> w = WedgeAtoms(0.1, molecule('CO'), pbc=True)
        >>> w.set_pbc([False,False,False]) # is an error in)
        """
        
        # protect the atoms.pbc from direct manipulation 
       
       
        # a True value, None would  harm here if _pbc_repr_int_not_bool is false
        pbc_radial_alway_True = 2 
        
        # for WedgeAtoms there is no third dimention describe periodicity
        #FIXME: impement numpy.nan or None insted of 0 for safety  
        pbc_not_an_atribute = 0    
        
        if isinstance(pbc, numpy_ndarray):
            #print 'In set_pbc -> ndarray inputed!'
            #print 'pbc[0]=',pbc[0]; print 'pbc[1]=',pbc[1]; print 'pbc[2]=',pbc[2]
            # FIXME: pbc_z attribute!
            pbc = pbc[1]
            #type(pbc) is 'numpy.bool_'now
            pbc = (pbc_radial_alway_True, pbc, pbc_not_an_atribute)
            
        elif isinstance(pbc, int):
            pbc = (pbc_radial_alway_True, pbc, pbc_not_an_atribute)
        else:
           raise TypeError('Method set_pbc(self, pbc): You may want to supply' 
                            + ' int type for pbc attribute '
                            + '(%s given)' % type(pbc))
        
        if self._pbc_repr_int_not_bool:
            self.pbc = numpy_array(pbc, int)
        else:
            self.pbc = numpy_array(pbc, bool)
   
#========== older functionality  ================================================
#       def set_pbc(self, pbc):
#        """
#        Use set_pbc_z() instead to set pbc along Z axis.
#        
#        For backward compatibility with ase.Atoms and 
#        should not be used unless you are not satisfied with set_pbc_z()
#        
#        (Example of not a proper usage:
#        >>> w = WedgeAtoms(0.1, molecule('CO'), pbc=True)
#        >>> w.set_pbc([False,False,False]) # is an error in)
#        """
#        
#        # None will harm here if _pbc_repr_int_not_bool is false
#        pbc_x_y_always_true = 2   
# 
#        if isinstance(pbc, (list, tuple)):
#            # FIXME: exception
#            #print 'in set_pbc'
#            #print 'pbc:', pbc
#            #print 'pbc.type:', type(pbc)
#            #print 'pbc.tolist:', pbc.tolist()
#            raise Exception
#        
#        if isinstance(pbc, int):
#            pbc = (pbc_x_y_always_true, pbc_x_y_always_true, pbc)
#        
#        
#        if type(pbc) == type(numpy_zeros(1)):
#           
#            #print 'in set_pbc -> ndarray inputed!'
#            #print 'pbc[0] = ', pbc[0]
#            #print 'pbc[1] = ', pbc[1]
#            #print 'pbc[2] = ', pbc[2]
#            # FIXME: z attribute!
#            pbc = pbc[2]
#            pbc = (pbc_x_y_always_true, pbc_x_y_always_true, pbc)
# otherwise it must be an instance of isinstance(pbc, (npy.ndarray,)):
#        # which is ok  
#        # or a bad user input 
#        
# #===============================================================================
# #        if isinstance(pbc, (npy.ndarray,))
# #            # FIXME: hardcoded data!
# #            pbc = True
# #            pbc = (pbc_x_y_always_true, pbc_x_y_always_true, pbc)
# #===============================================================================
# 
#        if self._pbc_repr_int_not_bool:
#            #self.pbc = npy.array(pbc, int)
#            self.pbc = numpy_array(pbc, int)
#        else:
#            self.pbc = numpy_array(pbc, bool)
#===============================================================================
            
       
    def get_pbc(self):
        """
        Use get_pbc_z() instead to set pbc along Z axis.
        
        For backward compatibility with ase.Atoms and 
        should not be used unless you are not satisfied with get_pbc_z()
        
        Get periodic boundary condition flags (in internal format)
        """
        # FIXME: generate warning instead of print
        #print 'Warning: %s.get_pbc() was used!' % self.__class__.__name__
        
        #return self.pbc.copy()
        
    
        return Atoms.get_pbc(self)
    
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
        
        #print self.pbc[0]
        #print self.pbc[1]
        #print bool(self.pbc[2])
        
        ##return bool(self.pbc[2])
        # periodicuty in z dirrection is decoded in second ([1]) varilbe  
        return bool(self.pbc[1])

#========= older functionality ===============================================
#    def get_pbc_z(self):
#        """
#        Get periodic boundary condition flag in z direction! 
#        (in internal format)
#        """
#        
#        #print self.pbc[0]
#        #print self.pbc[1]
#        #print bool(self.pbc[2])
#        
#        return bool(self.pbc[2])
#==============================================================================
        
    # DOC: this
    pbc_z = property(fget=get_pbc_z, fset=set_pbc_z,
                     doc='periodic boundary condition flag in z' 
                          + 'direction. Uses get/set methods')
    #pbc = property(fget=get_pbc_z, fset=set_pbc_z, doc='New pbc attribute!')
     
    def copy(self):
        # Why Atoms, not to use __init__ instead?!
        """
        Return a copy.
        """
        # TODOed: not tested -> seam to work 
        
        wedgeatoms = WedgeAtoms(angle=self.angle,
                                pbc_z=self.get_pbc_z(),
           # deep copy of all arrays, constraints, adsorbate_info...
           # must be done Atoms.__init__
                                symbols=Atoms.copy(self))
        
        return wedgeatoms

    
    def __repr__(self):
        """
        Return informal object representation 
        (see Comments (1.) in the class docstring )
        
        Note, accurate behavior relies ase.Atoms.__repr__(self)
        return, in particular,
        
        1. it starts from 'Atoms(' and 2. it contains pattern 
        'pbc=[ <zero or more of any character >]'     
        """
        # FIXME: accurate behavior relies ase.Atoms.__repr__(self) return  
        
        s = Atoms.__repr__(self)
        
        s = s[: - 1] + ', number_of_cells=%s)' % self.number_of_cells.tolist()
        
        import re
        pattern = r'pbc=\[.*\]' # pbc=[ <zero or more of any character >]
        if not (s.startswith('Atoms(') and re.search(pattern, s)):
            # FIXME: generate warning 
            s += '*Warning: this might be not accurate! If ' 
            s += 'ase.Atoms.__repr__() changes,'
            s += self.__class__.__name__
            s += '.__repr__() might need changes' 
        # add class name in front
        # insert  'angle = ...'  after 'Atoms(' at the begining  
        s = self.__class__.__name__ + '(angle = %s, ' % repr(self.angle) + s[6:]
        
        self._pbc_output_full = True
        if not self._pbc_output_full:
            # insert 'pbc_z= ...'
            s = re.sub(pattern, 'pbc_z=' + repr(self.get_pbc_z()), s)
        return s
    
    def __eq__(self, other):
        """
        Check for identity of two wedge atoms  objects.

        Identity means: same wedge angles, positions, atomic 
        numbers, unit cell and  periodic boundary conditions.
        """        
        
        return (Atoms.__eq__(self, other) and 
                self.angle == other.angle) 
        
        
    
    def _get_angle(self):
        """
        Return Wedge angle in radians
        """
        return self._angle
    
    def _set_angle(self, angle):
        # FIXME: what if angle/2_pi is not int?!
        # within tolerance
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
                raise Exception     
        else:
            raise ValueError('angle is not in required range (%f,%f]. Given %f)'
                                 % (self._min_angle, self._max_angle, angle))
        
    # Should it be a property ?!
    # with __slot__ there is not much point
    angle = property(fget=_get_angle, fset=_set_angle,
                     doc='Wedge angle in radians')
        
    
    def _is_wedge_atoms_instance(self):
        """
        Return True if self is WedgeAtoms instance 
        """
        return self._is_wedge_atoms_instance
        #return self._wedge_atoms_flag
    
    
    # FIXEME: name is_wedge_atoms
    # FIXME: What? _is_wedge_atoms_instance not a hidden argument!
    # isn't function not enough?
    is_wedge_atoms_instance = property(fget=_is_wedge_atoms_instance,
                                       doc='True if self is WedgeAtoms' 
                                       + 'instance (or subclass)')
    
    
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
    
    
    def get_max_radial_image_number(self):
        """
        Maximum nmber of radial images/cells in 
        """
    
        return int((2 * pi) / self.angle)
    

    def get_copies (self, to=None, from_=0, center_atoms=False):
        #FIXME: to -> (to -1)!!
        # FIXME: Return type Atoms?!
        # FIXME: name: get_rotated_copies
        # DOC:   (from_+1) is an id of first image to copied 
        """
        Return Atoms object containig all images of self by defeault, 
        i.e. if 'to' is None.
        
        Otherwise, if 'from_' is zero, return number of images
        specisiled by 'to'.  
        
        Whether 'from_' is zero or not, return images statring from 
        an image with id 'from_' to image with id ('to'+1), 
        so that ('to'-'from_') is a number of images returned.  
        
        See get_copy() for more info on image definition and indexing.
        
        @param to: int
            (to-1) an id of last image to be copied
        @param from_: int 
            from_ is an id of first image to copied
            from_ should be smaller then to
        @param center_atoms: bool 
            wheather to center return atoms w.r.t unit cell
            (translate atoms so that center mass 
            is in the center of the unit cell)
            
        """
        
        # FIXME: Implement feature 
        if center_atoms:
            raise NotImplementedError('Please, use center_atoms=Flase')
        
      
        id_range = self.get_max_radial_image_number()
                
        if to is None:
            to = id_range
               
        if (not isinstance(to, int)) or (not isinstance(from_, int)):
            raise TypeError('Both to and from_ attributes should be int')
            
        if (from_ >= to):
            raise ValueError('from_ (2nd attribute) should be strictly'
                             + ' smaller then to (1st attribute).' 
                             + ' Given: from_=%d, to=%d.' % (to, from_))
        
        if (to - from_) > id_range:
            raise ValueError('For current wedge angle (%f),' % (self.angle) 
                             + ' maximum image number could be % d,' 
                             % (id_range) 
                             + ' while to-from_ = %d.' % (to - from_))
        a = Atoms() 
        for id in range(from_, to):
            a += self.get_copy(id)
        return a 
    
    
    
    
    
    def get_copy (self, id=0):
        #FIXME: add all new attributes!!
        # DOC: image -> radial/rotational image/ OR cell 
        # FIXME: name: get_rotated_copy
        # FIXME: Return type Atoms?! 
        #        insted return WedgeAngel, w. angle_from?
        """
        Return Atoms object containig an images of self with given id.
        
        An image (id ==1) is a self rotated counterclockwise 
        by the angel 'self.angle' around 'z' axis.
        
        id=0 is equivalent to self (within return type)   
        id=-1 is equivalent to self rotated clockwise 
            by the angel 'self.angle'
        id=1 is equivalent to self rotated 
            counterclockwise by the angel 'self.angle'
        id=2 is equivalent to self rotated 
            counterclockwise twice by the angel 'self.angle'
        ...
               
        @param id: int, (Defalut: 0)
            image id, e.g.1, +2, -1.
            Can also be float, eg id=0.5 rotates self 
            by half of the angel 
            (Use float values with caution!)
        """
        if not isinstance(id, (int, float)):
            raise TypeError('id (the only attribute) should be int or float')
        
        atoms = Atoms(self)
        rotate_angle = id * self.angle
        atoms.rotate('z', rotate_angle)
        return atoms
                 
        #return self._get_radial_image_atoms(Atoms(self), id)

#===============================================================================
#    def _get_radial_image_atoms(self, atoms, radial_image_id):
#        # DOC:
#        """
#        As get_copy, but applies to atoms 
#        """
#                
#        rotatition_angle = radial_image_id * self.angle
#        atoms.rotate('z', rotatition_angle) 
#        return atoms  
#===============================================================================
    
    
#===============================================================================
#    def _get_vector_from_radial_image_oldest(self, atom_index, radial_image_id):
#        # FIXME: more effective implementaion that makes rotate  
#        """
#        Return coordinates of atom atom_index in radial image 
#        radial_image_id 
#        
#        Not effective for image_id = 0 and should not be used for it 
#        """
#        # FIXME: Atoms([self[atom_index]])
#        atom_in_image = self._get_radial_image_atoms(
#                    Atoms([self[atom_index]]), radial_image_id)
#        
#        #print 'atom_in_image.positions', atom_in_image.positions 
#        
#        return atom_in_image.positions[0]
#    
#===============================================================================
    
#===============================================================================
#    def _get_vector_from_radial_image_older(self, atom_index, radial_image_id):  
#        """
#        ...
#        
#        # with resepect to origin!!
#         
#        """
#        # with resepect to origin!!
#        
#        position = self[atom_index].position
#        x_old = position[0]
#        y_old = position[1]
#        z_old = position[2]
#        
#        ro = numpy_linalg_norm(numpy_array([x_old, y_old]))
#        #ro2 = sqrt(x_old ** 2. + y_old ** 2.)
#        #if ro != ro2:
#        #    print 'error!!'
#        #    print ro
#        #    print ro2
#        #else:
#        #    print 'Ok!'
#        cos_old = x_old / ro
#        sin_old = y_old / ro
#        
#        rotatition_angle = radial_image_id * self.angle
#        cos_rot = cos(rotatition_angle)
#        sin_rot = sin(rotatition_angle)
#        
#        cos_new = cos_old * cos_rot - sin_old * sin_rot
#        sin_new = sin_old * cos_rot + cos_old * sin_rot
#    
#        x_new = ro * cos_new
#        y_new = ro * sin_new
#        
#        #FIXME: typecheck
#        dtype = self[atom_index].position.dtype
#        new_positions = numpy_array([x_new, y_new, z_old], 
#                                    dtype = dtype) 
#        
#        #print "2:", new_positions
#        #print "1:", self._get_vector_from_radial_image(atom_index, radial_image_id)
#        
#        return new_positions
#===============================================================================
    
#===============================================================================
#    def _get_vector_from_z_image(self, atom_index, z_image_id):
#        """
#        Critically depends on correctn unit cell settings, 
#        on dz in particular
#        """
#        position = self[atom_index].position
#        
#        translate_vector = z_image_id * self.get_cell()[2]
#        
#        moved_position = position + translate_vector
#        #print 'In _get_vector_from_z_image'
#        #print 'Original_vector = ', position
#        #print 'Moved_vector= ', moved_position
#               
#        return moved_position
#===============================================================================

    
    
    def _get_vector_from_radial_image(self, position, radial_image_id):  
        # FIXME: dirrect implemenationmn of rotate 
        """
        Rorate the position vector in xy-plane 
        for rotatition_angle = radial_image_id * self.angle
        
        # with resepect to origin!!
        """
        
        x_old = position[0]
        y_old = position[1]
        z_old = position[2]
        
        ro = numpy_linalg_norm(numpy_array([x_old, y_old]))
        
        cos_old = x_old / ro
        sin_old = y_old / ro
        
        rotatition_angle = radial_image_id * self.angle
        cos_rot = cos(rotatition_angle)
        sin_rot = sin(rotatition_angle)
        
        cos_new = cos_old * cos_rot - sin_old * sin_rot
        sin_new = sin_old * cos_rot + cos_old * sin_rot
    
        x_new = ro * cos_new
        y_new = ro * sin_new
        
        new_positions = numpy_array([x_new, y_new, z_old],
                                    dtype=numpy_float64) 
               
        return new_positions
    

        
    def get_vector(self, i, j, radial_mic=True, z_mic=False):
        # FIXME: MIC = Minimum Image Convention 
        # or minimum image criterion?!
        # http://www.fisica.uniud.it/~ercolessi/md/md/node18.html#SECTION00421000000000000000
        """
        Return vector between two atoms (R_{ij} = R_j - R_i).  
        
        Default is to use the Minimum Image Convention (radial_mic=True) 
        
        MIC works withing WedgeAtoms images only!
        
        
        See get_copy() for more info on image definition and indexing.   
        """
            
        verbouse = False
        
        if z_mic:
            raise NotImplementedError('z_mic=True is not implemented yet')
        
        R_i = self.arrays['positions'][i]
        R_j = self.arrays['positions'][j]
        
        vector = numpy_zeros((3, 3))
        vector[0] = R_j - R_i
        
        min_index = 0
        
        if radial_mic: 
            
            # vector[1] between images +1 and 0
            vector[1] = self._get_vector_from_radial_image(R_j, + 1) - R_i
            # vector[2] between images -1 and 0
            vector[2] = self._get_vector_from_radial_image(R_j, - 1) - R_i 
            
            
            _tolerance = 0.00001
            # serch for min norm of vector[0], vector[1], vector[2]            
            for index in range(3):
                norm = numpy_linalg_norm(vector[index])
                min_norm = numpy_linalg_norm(vector[min_index])
                 
                if verbouse:
                    #print "V[%d]= " % index, vector[index]
                    print "Norm[%d] = " % index, norm
                
                # FIXME: no tolerance needed?
                if (min_norm - norm > _tolerance): 
                    min_index = index
            if verbouse:
                print "min_index = ", min_index
                #print "min norm = ", numpy_linalg_norm(vector[min_index])
                print "min vector = ", vector[min_index]
    
        return vector[min_index]
    
    def get_distance(self, a0, a1, mic=False):
        """
        Return distance between two atoms.
                
        Default is to *not* use the Minimum Image Convention, 
        mic=False (for ase backward compatibility)
        
        mic=True refers to *radial only* Minimum Image Convention, 
        which works for WedgeAtoms only! 
        """
        vector = self.get_vector(a0, a1, radial_mic=mic)
        
        return numpy_linalg_norm(vector)
    
    
    # DOC:
    def get_number_of_cells(self):
        return self._number_of_cells
    
    
    # DOC:
    def set_number_of_cells(self, number_of_cells=None):
        """
        is None then according to pbc_z and max radial range   
        
        @param number_of_cells:
        
        set_number_of_cells((10, numpy.Inf, 1))
        set_number_of_cells((10, numpy.Inf))
        
        
        if pbc_z = Inf
        """
        
        if number_of_cells is None:
            if self.get_pbc_z():
                z_default = numpy_Inf
            else: 
                z_default = 1
            number_of_cells = (self.get_max_radial_image_number(),
                                   z_default)
        #FIXME: is this good?!
        if isinstance(number_of_cells, int):
            #number_of_cells = (number_of_cells, number_of_cells)
            raise TypeError('Tuple or list required')
        
        if isinstance(number_of_cells, (list, tuple)):
            if len(number_of_cells) == 2:
                number_of_cells = (number_of_cells[0],
                                   number_of_cells[1],
                                   1)
            elif len(number_of_cells) != 3:
                raise ValueError('number_of_cells tuple/list should' + 
                                 ' containt 2 of 3 items. ' 
                                 + '(%d supplied)' % len(number_of_cells)) 
        
        for value in number_of_cells:
            if value <= 0:
                raise ValueError('Each number_of_cells values should '
                                 + 'be greater thet than 0. '
                                 + ' (Supplied %d)' % (value))

        
        #print number_of_cells[0]
        #print self.get_max_radial_image_number()
        if number_of_cells[0] > self.get_max_radial_image_number():
            raise ValueError('For current wedge angle (%f),' % (self.angle) 
                             + ' maximum number of cells in radial'
                             + ' dirrection could be '
                        + '%d' % (self.get_max_radial_image_number())
                             + ' (supplied %d).' % (number_of_cells[0]))
        
        
        if number_of_cells[2] != 1:
            raise ValueError('Third value (if given) should = 1, '  
                             + ' since for WedgeAtoms there is no third' 
                             + ' dimention describe periodicity')
        
        
        # FIXME: add Inf
        self._number_of_cells = numpy_array(number_of_cells)
    
        
    # FIXME: add to __rerp__ 
    # FIXME: add to __init__
    number_of_cells = property(fget=get_number_of_cells,
                               fset=set_number_of_cells,
                               doc='Shotcut for set_number_of_cells()'  
                                    + 'and get_number_of_cells')
    
    
    # older symmetry_operation
        #----------------------- def symmetry_operation(self, atom_index, image_id):
        #------- #TODO: cahnge everywhere image_id to cell_id or image_id or vv.
        #---------------------------------------- #TODO: cell or image of cell?!
        #------------------------------------------------------------------- """
        #------------------------------ Return vector characterizing position of
        #--------------------------------- atom atom_index in cell with cell_id.
        #------------------------------------------------------------------------------ 
        #---------------------------------------------------- @param atom_index:
        #------------------------------------------------------ @param image_id:
            #---------------------------------- (radial_image_id, z_image_id) or
            #----------------------------------- (radial_image_id, z_image_id,0)
            #------------------------------------------------------------------------------ 
        #-------------------------------------------------------------- Example:
        #------------------------------------- symmetry_operation(5, (1, -5, 0))
        #---------------------------------------- symmetry_operation(5, (1, -5))
        #------------------------------------------------------------------- """
        #------------------------------------------------------------------------------ 
        #------------------------------------------------------------------------------ 
        #----------------------------------------- if isinstance(image_id, int):
            #------------------------- raise TypeError('Tuple or list required')
            #------------------------------------------------------------------------------ 
        #------------------------------- if isinstance(image_id, (list, tuple)):
            #-------------------------------------------- if len(image_id) == 2:
                #-------------------------------------- image_id = (image_id[0],
                            #-------------------------------------- image_id[1],
                            #------------------------------------------------ 0)
            #------------------------------------------ elif len(image_id) != 3:
                #--------------- raise ValueError('image_id tuple/list should' +
                                 #------------------- ' containt 2 of 3 items. '
                                 #---------- + '(%d supplied).' % len(image_id))
                                 #------------------------------------------------------------------------------ 
        #--------------------------------------------------- if image_id[0] < 0:
            #------------------ raise ValueError('image_id  in radial dirrecton'
                             #-- + '(first attribute) could not be less then 0 '
                             #------------- + ' (Supplied %d).' % (image_id[0]))
                             #------------------------------------------------------------------------------ 
        #---------- # FIXME: add chec so that they are in number_of_cells range!
        #------------------------------------------------------------------------------ 
        #------------------------- # max image_id == self.number_of_cells()[0]-1
        #---------------------------------- #print self.get_number_of_cells()[0]
        #---------------------- if image_id[0] >= self.get_number_of_cells()[0]:
            #------------- raise ValueError('For current number_of_cells value,'
                             #------------------ + ' maximum image id in radial'
                             #------ + ' dirrection (first attribute) could be '
                             #----- + '%d' % (self.get_number_of_cells()[0] - 1)
                             #------------- + ' (supplied %d).' % (image_id[0]))
                             #------------------------------------------------------------------------------ 
        #-------------------------------------------------- if image_id[2] != 0:
            #------------ raise ValueError('Third value (if given) should = 0, '
                             #------ + ' since for WedgeAtoms there is no third'
                             #----------- + ' dimention describing periodicity')
                             #------------------------------------------------------------------------------ 
        #------------------------------------------ # FIXME: range check for z?!
        #------------------------------------------------------------------------------ 
        #----------------------------------------- radial_image_id = image_id[0]
        #---------------------------------------------- z_image_id = image_id[1]
        #------------------------------------------------------------------------------ 
        #------------------------------------------------------------------------------ 
        #---------------------------------- position = self[atom_index].position
        #------------------------------------------------------------------------------ 
        #----- #position_z_zero = self._get_vector_from_radial_image(atom_index,
        #--------- #                                            radial_image_id)
        #------------------------------------------------------------------------------ 
        #-------- position_z_zero = self._get_vector_from_radial_image(position,
                                                    #---------- radial_image_id)
                                                    #------------------------------------------------------------------------------ 
        #-------------------- translate_vector = z_image_id * self.get_cell()[2]
        #--------------------- new_position = position_z_zero + translate_vector
        #------------------------------------------------------------------------------ 
        #--------------------------------- #print 'Original_vector = ', position
        #---------------------------------- #print 'new_position:', new_position
        #------------------------------------------------------------------------------ 
        #--------------------------------------------------- return new_position
    
    #FIXME: chahe name to transforms 
    def symmetry_operation(self, atom, image_id):
        #FIXME: implement using the rotation matrix multiplications?! 
        #TODO: cahnge everywhere image_id to cell_id or image_id or vv. 
        #TODO: cell or image of cell?!
        """
        Return vector characterizing position of 
        atom atom_index in cell with cell_id.  
        
        atom is tuple -> position_vector
        atom is int -> atom_index
        
        @param atom_index:
        @param image_id:
            (radial_image_id, z_image_id) or
            (radial_image_id, z_image_id,0)
        
        Example:
        symmetry_operation(5, (1, -5, 0))
        symmetry_operation(5, (1, -5))
        """
        
        
        if isinstance(image_id, int):
            raise TypeError('Tuple or list required')
        
        if isinstance(image_id, (list, tuple)):
            if len(image_id) == 2:
                image_id = (image_id[0],
                            image_id[1],
                            0)
            elif len(image_id) != 3:
                raise ValueError('image_id tuple/list should' + 
                                 ' containt 2 of 3 items. ' 
                                 + '(%d supplied).' % len(image_id)) 
        
        if image_id[0] < 0:
            raise ValueError('image_id  in radial dirrecton'
                             + '(first attribute) could not be less then 0 '
                             + ' (Supplied %d).' % (image_id[0]))
        
        # FIXME: add chec so that they are in number_of_cells range!
             
        # max image_id == self.number_of_cells()[0]-1
        #print self.get_number_of_cells()[0]
        if image_id[0] >= self.get_number_of_cells()[0]:
            raise ValueError('For current number_of_cells value,'
                             + ' maximum image id in radial'
                             + ' dirrection (first attribute) could be '
                             + '%d' % (self.get_number_of_cells()[0] - 1)
                             + ' (supplied %d).' % (image_id[0]))
        
        if image_id[2] != 0:
            raise ValueError('Third value (if given) should = 0, '  
                             + ' since for WedgeAtoms there is no third' 
                             + ' dimention describing periodicity')
            
        # FIXME: range check for z?!
       
        radial_image_id = image_id[0]
        z_image_id = image_id[1]
         
        
        if isinstance(atom, int):
            position = self[atom].position
        else: 
            position = atom
        #elif (isinstance(atom, tuple) and len(atom) == 3):
        #    position = atom
        #else:
        #    raise TypeError("'atom' attribute should be either" +  
        #        "atom index or atom position (in 3d space)")
        
        
        #position_z_zero = self._get_vector_from_radial_image(atom_index,
        #                                            radial_image_id) 
        
        position_z_zero = self._get_vector_from_radial_image(position,
                                                    radial_image_id) 
        
        translate_vector = z_image_id * self.get_cell()[2]
        new_position = position_z_zero + translate_vector
        
        #print 'Original_vector = ', position
        #print 'new_position:', new_position
        
        return new_position    
    
    def transform_der(self, image_id):
        """
        Derevative of symmetry transfomation =  
        = Rotation matix, if transfomiont contains rotations  
        
        @param image_id:
        
        """
        # FIXME: add type checks for image_id 
        
        a = image_id[0]*self.angle
        #print image_id[0]
        #print a  
        
        #rotational matrix around z axis for active rotations
        #FIXME: from scipy?!  
        T = numpy_array([[cos(a), -sin(a), 0],
                         [sin(a),  cos(a), 0],
                         [     0,       0, 1]
                         ])
        return T
    
    
#===============================================================================
#    def get_cell(self):
#        """
#        for testing of drawing...only!
#        
#        Get the three unit cell vectors as a 3x3 ndarray.
#        """
#        a = self.cell.copy()
#        
# #        print;  print 'in get_cell' 
# #        print a
# #        
# #        a = numpy_array([[10., numpy_Inf,0.],
# #                   [7.07, 7.07106781,0.],
# #                   [0., 0., 0.]], 
# #                   dtype = numpy_float64)
# #        
#        return a
#===============================================================================


    def repeat(self, rep):
        """
        Works in z direction only
        
        
        rep = (radial_rep, z_rep, 1)
        rep = (radial_rep, z_rep)
        """
        
        print 'rep:', rep
        
        if isinstance(rep, int):
            #rep = (rep, rep)
            raise TypeError('Tuple or list required')
        
        if isinstance(rep, (list, tuple)):
            if len(rep) == 2:
                rep = (rep[0],
                       rep[1],
                       1)
            elif len(rep) != 3:
                raise ValueError('rep tuple/list should' + 
                                 ' containt 2 of 3 items. ' 
                                 + '(%d supplied)' % len(rep)) 
        
        for value in rep:
            if value <= 0:
                raise ValueError('Each rep values should '
                                 + 'be greater thet than 0. '
                                 + ' (Supplied %d)' % (value))

        
        #print rep[0]
        #print self.get_max_radial_rep()
        if rep[0] > self.get_max_radial_image_number():
            raise ValueError('For current wedge angle (%f),' % (self.angle) 
                             + ' maximum number of cells in radial'
                             + ' dirrection could be '
                        + '%d' % (self.get_max_radial_image_number())
                             + ' (supplied %d).' % (rep[0]))
        
        
        if rep[2] != 1:
            raise ValueError('Third value (if given) should = 1, '  
                             + ' since for WedgeAtoms there is no third' 
                             + ' dimention describing periodicity')
        
        
        #FIXME: regenerate WedgeAtoms!
        #FIXME: add all propetries
        
        wa = WedgeAtoms(self.angle,
                        Atoms.repeat((self.get_copies(rep[0])), (1,1,rep[1])))
        
        return wa  

    __mul__ = repeat
        
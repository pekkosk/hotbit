# Module under development!
 
from ase import Atoms as AtomsOriginal
from wedgeatoms import WedgeAtoms

class Atoms(WedgeAtoms):
    # TODO: not tested!
    """ 
    (Generalized) Atoms Object (extended ase.Atoms with WedgeAtoms)
    """
    
    def __init__(self,
                 # New arguments
                 angle=None,
                 pbc_z=None,
                 # old argument
                 pbc=None,
                 # TODO: cell might need to be changed
                 cell=None,
                 # Old arguments, same as in ase.Atoms 
                 symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 constraint=None,
                 calculator=None):
        
        """
        If angle is set not given, behaves as ase.Atoms
        Otherwise, behaves as WedgeAtoms.
        
        ''pbc'' is primary parameter. If both ''pbc'' and ''pbc_z'' given, 
        the value of ''pbc_z'' is discarded. Either of them carries 
        the meaning of ''pbc'' from ase.Atoms.__init__() or of ''pbc_z'' 
        from WedgeAtoms__init__(). 
        (specific behavior is dependent on ''angle'' value).    
        """
        
        if pbc is None: 
            pbc = pbc_z
        # otherwise value of pbc_z is discarded 
        
        if angle is None:
            # behaves as typical ase.Atoms
            AtomsOriginal.__init__(self, pbc=pbc, cell=cell,
                 symbols=symbols,
                 positions=positions, numbers=numbers,
                 tags=tags, momenta=momenta, masses=masses,
                 magmoms=magmoms, charges=charges,
                 scaled_positions=scaled_positions,
                 constraint=constraint,
                 calculator=calculator)
        else: 
            # behaves as WedgeAtoms
            WedgeAtoms.__init__(self, angle=angle,
                 pbc_z=pbc,
                 # TODO: cell might need to be changed
                 cell=cell,
                 symbols=symbols,
                 positions=positions, numbers=numbers,
                 tags=tags, momenta=momenta, masses=masses,
                 magmoms=magmoms, charges=charges,
                 scaled_positions=scaled_positions,
                 constraint=constraint,
                 calculator=calculator)
        
    
    

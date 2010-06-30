# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

import numpy as np
from ase.units import Hartree, Bohr
from weakref import proxy

def setup_vdw(calc):
    if calc.get('vdw'):
        raise NotImplementedError('van der Waals interactions are not yet implemented.')
    
    elms = len(calc.el.present)
    for i,s1 in enumerate(calc.el.present):
        for s2 in calc.el.present[i:]:      
            #
            # here vdw is simply the interaction
            # between elements s1 and s2
            #                 
            calc.pp.add_pair_potential(s1,s2,vdw,eVA=True)
    

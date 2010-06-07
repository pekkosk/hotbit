from hotbit import Hotbit
from ase import *
from hotbit.test.misc import default_param
from hotbit.test.misc import molecule
from numpy import *
from box.md import microcanonical, energy_conservation
from ase.units import fs
import sys

from box.fd_forces import check_forces

systems=[
    'AuC',
    'H2COH'
    ]

#default_param['convergence']=1E-5
#default_param['Anderson_memory']=3
default_param['width']=0.1

for charge_density in [ None, 'Gaussian', 'Slater' ]:
    for system in systems:
        if charge_density is None:
            print '    ... forces for %s, no SCC' %system
            calc=Hotbit(verbose=True,SCC=False,txt='forces.cal',**default_param)
        else:
            print '    ... forces for %s, SCC, charge density = %s' % \
                ( system, charge_density )
            calc=Hotbit(verbose=True,SCC=True,charge_density=charge_density,
                        txt='forces.cal',**default_param)
        atoms=molecule(system)
        #view(atoms)
        atoms.center(vacuum=5)
        atoms[0].x+=0.1
        atoms=Atoms(atoms)
        atoms.set_calculator(calc)
        sys.stdout.flush()
        rec, de, du=energy_conservation(atoms,dt=0.2,steps=100)
        assert de/du<0.011
        
    





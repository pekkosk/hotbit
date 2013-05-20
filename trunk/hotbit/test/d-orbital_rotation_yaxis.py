from ase import *
from hotbit import *
from ase import io

import numpy as np
from box.mix import phival
from math import sin,cos 
from weakref import proxy
import warnings


M=7

atoms = Atoms('Au2',[(5,0,0),(5,2.5,0.3)],container='Wedge')


atoms.set_container(M=M,height=5.0)
calc=Hotbit(SCC=False,txt='-',kpts=(M,1,1))
atoms.set_calculator(calc)
e1 = atoms.get_potential_energy()


atoms = Atoms('Au2',[(0,0,5),(2.5,0.3,5)],container='WedgeYAxis')
atoms.set_container(M=M,height=5.0)
calc=Hotbit(SCC=False,txt='-',kpts=(M,1,1))
atoms.set_calculator(calc)
e2 = atoms.get_potential_energy()

assert np.abs(e1-e2)<1E-13

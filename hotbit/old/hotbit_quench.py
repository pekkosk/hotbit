#!/usr/bin/env python
"""
    Hotbit structure optimizer script.
    
    Usage:
        hotbit_quench.py file.xyz
"""
from hotbit import Calculator
from ase import *
from hotbit.misc import quench
from box import mix
import sys

try:
    file=sys.argv[1]
except:
    print __doc__
    sys.exit()
    
name=file.split('.')[0]
atoms=read(file)
calc=Calculator()
atoms.set_calculator(calc)
quench(atoms,name=name,fmax=0.05,calc=calc)

calc.finalize()

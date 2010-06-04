from ase import *
from hotbit import Hotbit
import numpy as nu
import os

dir = os.environ['HOTBIT_DIR']+'/param/'

atoms = Atoms('CH2',((0,0,0),(1,0,0),(-1,0,0)))
atoms.center(vacuum=5)
calc = Hotbit(charge=-2, txt='-',
    elements={'H':dir + 'H.elm', 'C':dir+'C.elm'},
    tables = {'HH':dir+'H_H.par', 'CC':dir+'C_C.par', 'CH':dir+'C_H.par'},
    SCC=False,
    width=0.04,
    mixer='pulay')
atoms.set_calculator(calc)
e1 = atoms.get_potential_energy()

calc.write()
calc2 = Hotbit(txt='-', filename='restart.hb')
atoms2 = calc2.get_atoms()
e2 = atoms2.get_potential_energy()

if abs(e1-e2) > 1e-7:
    raise RuntimeError("The original and loaded calculator doesn't give same results!")

atoms = Atoms('CH2',((0,0,0),(1,0,0),(-1,0,0)),pbc=[1,1,1])
atoms.center(vacuum=5)
calc = Hotbit(charge=0, txt='-',
    elements={'H':dir + 'H.elm', 'C':dir+'C.elm'},
    tables = {'HH':dir+'H_H.par', 'CC':dir+'C_C.par', 'CH':dir+'C_H.par'},
    SCC=False,
    width=0.04,
    mixer='pulay')
atoms.set_calculator(calc)
e1 = atoms.get_potential_energy()

calc.write()
calc2 = Hotbit(txt='-', filename='restart.hb')
atoms2 = calc2.get_atoms()
e2 = atoms2.get_potential_energy()

if abs(e1-e2) > 1e-7:
    raise RuntimeError("The original and loaded calculator doesn't give same results!")
os.system('rm -f restart.hb')

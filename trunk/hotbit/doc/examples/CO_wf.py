from ase import *
from hotbit import *
from hotbit import fixpar
import os

calc=Hotbit(SCC=True,width=0.05,txt='test.cal',parameters=fixpar)
atoms=Atoms('CO',positions=[(0,0,0),(1.13,0,0)],pbc=False)
atoms.center(vacuum=3) 
atoms.set_calculator(calc)
atoms.get_potential_energy()

wf=calc.get_grid_wf(0,spacing=0.5)
write('wf0.cube',atoms,data=wf)

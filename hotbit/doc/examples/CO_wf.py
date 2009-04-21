from ase import *
from hotbit import Calculator
from hotbit.analysis import WaveFunctions

calc=Calculator(SCC=True,width=0.05,txt='test.cal')
atoms=Atoms('CO',positions=[(0,0,0),(1.13,0,0)],pbc=False)
atoms.center(vacuum=3) 
atoms.set_calculator(calc)
atoms.get_potential_energy()

wf=WaveFunctions(atoms,dr=0.5) #0.5 Angstrom grid
wf0=wf.get_wf(0)
write('wf0.cube',atoms,data=wf0)

 # calculate CO dimer
from ase import *
from hotbit import Hotbit

for SCC in [False,True]:
    calc=Hotbit(SCC=SCC,width=0.05,txt='test.cal')
    atoms=Atoms('CO',positions=[(0,0,0),(1.13,0,0)],cell=(10,10,10),pbc=False)
    atoms.center(vacuum=10) #not necessary
    atoms.set_calculator(calc)
    print(atoms.get_potential_energy())

 # calculate Au dimer
from ase import *
from hotbit import Hotbit

calc=Hotbit(SCC=True,width=0.05,txt='test.cal')
Au2=Atoms('Au2',positions=[(0,0,0),(2.6,0,0)],cell=(10,10,10),pbc=False)
Au2.center(vacuum=10)
Au2.set_calculator(calc)
print(Au2.get_potential_energy())
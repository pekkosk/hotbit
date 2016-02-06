from ase import *
from hotbit import *
from ase.data.molecules import molecule

atoms = molecule('C6H6')
calc = Hotbit(SCC=True,width=0.05,txt='benzene.cal')
atoms.set_calculator(calc)
print atoms.get_potential_energy()
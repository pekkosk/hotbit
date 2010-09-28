from ase import *
from ase import Atoms as ase_Atoms
from hotbit import *
from hotbit.atoms import Atoms
from box.md import check_energy_conservation
from hotbit.test.misc import default_param

SCC=True
cut=3.0
atoms = Atoms('CH',[(1.42,0,0),(2.42,0,0)],container='Wedge')
atoms.set_container(M=6,height=10)

calc = Hotbit(SCC=SCC,txt='tmp.cal',kpts=(6,1,1),gamma_cut=cut,**default_param)
atoms.set_calculator(calc)
e1 = atoms.get_potential_energy()

 
aux = atoms.extended_copy((3,1,1))
atoms2 = Atoms(container='Wedge')
atoms2 += aux[0]
atoms2 += aux[-1]
atoms2.set_container(M=6,height=10)

#view(atoms2)
calc = Hotbit(SCC=SCC,txt='tmp.cal',kpts=(6,1,1),gamma_cut=cut,**default_param)
atoms2.set_calculator(calc)
e2 = atoms.get_potential_energy()

assert abs(e1-e2)<1E-5


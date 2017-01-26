from ase import *
from ase import Atoms as ase_Atoms
from ase.units import fs
from hotbit import *
from hotbit.atoms import Atoms
from box.md import check_energy_conservation
from hotbit.test.misc import default_param


# check that C1H1-presentation of C6H6 goes right
SCC=True
cut=3.0
atoms = Atoms('CH',[(1.42,0,0),(2.0,1.0,0.2)],container='Wedge')
atoms.set_container(M=6,height=10)

calc = Hotbit(SCC=SCC,txt='tmp.cal',kpts=(6,1,1),gamma_cut=cut,**default_param)
atoms.set_calculator(calc)
e1 = atoms.get_potential_energy()

atoms6 = ase_Atoms(pbc=False)
atoms6 += atoms.extended_copy([(i-2,0,0) for i in range(6)])
#view(atoms)
calc = Hotbit(SCC=SCC,txt='tmp.cal',gamma_cut=cut,**default_param)
atoms6.set_calculator(calc)
e6 = atoms6.get_potential_energy()

assert abs(6*e1-e6)<1E-5


#
# energy conservation
#
atoms = Atoms('CH',[(1.42,0,0),(2.0,0.5,0.3)],container='Wedge')
atoms.set_container(M=6,height=10)
calc = Hotbit(SCC=SCC,txt='tmp.cal',kpts=(6,1,1),gamma_cut=cut,**default_param)
atoms.set_calculator(calc)

atoms.rattle(0.1)
assert check_energy_conservation(atoms,dt=0.2*fs,steps=30,tol=0.01,plot=False)


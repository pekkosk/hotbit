from ase import *
from ase.units import fs
from ase.io import Trajectory
from hotbit import *
from box.md import check_energy_conservation
from ase import io

M=7
atoms = Atoms('Au2',[(5,0,0),(5,2.5,0.3)],container='Wedge')
atoms.set_container(M=M,height=5.0)
calc=Hotbit(SCC=False,txt='-',kpts=(M,1,1))
atoms.set_calculator(calc)

e1 = atoms.get_potential_energy()

whole = atoms.extended_copy((M,1,1))
calc=Hotbit(SCC=False,txt='-')
whole.set_calculator(calc)
e2 = whole.get_potential_energy()

assert abs(M*e1-e2)<1E-10


assert check_energy_conservation(atoms,dt=2.5*fs,steps=50,tol=0.01,plot=False)

#dyn = VelocityVerlet(atoms,2.5*fs)
#traj = Trajectory('koe.traj','w',atoms)
#dyn.attach(traj.write)
#dyn.run(100)

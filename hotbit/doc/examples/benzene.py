from ase import *
from ase.io import Trajectory
from hotbit import *
from ase.data.molecules import molecule

atoms=molecule('C6H6')
atoms.set_pbc(False)
atoms.center(vacuum=3)
calc=Hotbit(SCC=True, txt='benzene.cal')
atoms.set_calculator(calc)

traj=Trajectory('quench.traj','w',atoms=atoms)
qn=QuasiNewton(atoms)
qn.attach(traj.write)
qn.run()

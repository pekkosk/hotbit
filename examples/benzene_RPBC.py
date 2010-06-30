from ase import *
from hotbit import *
from ase.data.molecules import molecule

C6H6 = molecule('C6H6')
calc = Hotbit(SCC=True,width=0.05,txt='benzene.cal')
C6H6.set_calculator(calc)
e6 = C6H6.get_potential_energy()


atoms = Atoms(container='Wedge')
atoms += C6H6[0]
atoms += C6H6[6]
atoms.set_container(angle=2*pi/6,height=5,pbcz=False)

calc = Hotbit(SCC=True,width=0.05,txt='benzene_RPBC.cal',kpts=(6,1,1))
atoms.set_calculator(calc)
e1 = atoms.get_potential_energy()
print e6/6,' should equal',e1

traj=PickleTrajectory('quench.traj','w',atoms=atoms)
qn=QuasiNewton(atoms)
qn.attach(traj.write)
qn.run()
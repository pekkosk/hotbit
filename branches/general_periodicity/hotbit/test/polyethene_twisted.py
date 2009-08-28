from ase import *
from hotbit import *
from numpy import *
from box.md import check_energy_conservation

d=1.4
atoms = Atoms('C4H4',[(0,0,0),(0.5,0,d),(0,0,2*d),(0.5,0,3*d),\
                      (-1,0,0),(1.5,0,d),(-1,0,2*d),(1.5,0,3*d)],\
                      container='Chiral')
atoms.set_container(angle=2*pi/30,height=4*d)
#view(atoms)

calc=Hotbit(SCC=True,txt='polyethene.cal',gamma_cut=3*d,kpts=(1,1,20))
atoms.set_calculator(calc)
q = FIRE(atoms,trajectory='polyethene.trj',logfile=None)
q.run(fmax=0.5)

assert check_energy_conservation(atoms,dt=0.5*fs,steps=50,tol=0.01,plot=False)
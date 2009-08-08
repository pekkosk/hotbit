from ase import *
from hotbit import *
from hotbit.atoms import Atoms
from box.systems import nanotube
from box.md import check_energy_conservation

nkpts=10

# energy of normal infinite (5,0) 
straight = nanotube('C',1.42,5,0)
straight.set_pbc((False,False,True))
calc = Hotbit(SCC=False,txt='chiral.cal',kpts=(1,1,nkpts))
straight.set_calculator(calc)
e1 = straight.get_potential_energy()
#view(straight)


# same thing, but calculate by twisting 2*pi/5 while translating
height = straight.get_cell()[2,2]
chiral = Atoms(container='Chiral')
chiral += straight
chiral.set_container(height=height,angle=2*pi/5)
calc = Hotbit(SCC=False,txt='chiral.cal',kpts=(1,1,nkpts))
chiral.set_calculator(calc)
e2 = chiral.get_potential_energy()
assert abs(e1-e2)<1E-6


# check the energy conservation for the chiral situation
chiral.rattle(0.1)
calc = Hotbit(SCC=False,txt='chiral.cal',kpts=(1,1,1))
chiral.set_calculator(calc)
conv = check_energy_conservation(chiral,dt=0.5*fs,steps=50,tol=0.01,plot=False)
assert conv
from ase import *
from hotbit import *
from hotbit.atoms import Atoms
from box.systems import nanotube
from box.md import check_energy_conservation

nkpts=10
SCC=True
cut=5.0
# check that SCC works for chiral atoms

# energy of normal infinite (5,0) 
straight = nanotube('C',1.42,5,0)
straight.set_pbc((False,False,True))
for i in range(len(straight)):
    if straight[i].z<1.0:
        r=sqrt(straight[i].x**2+straight[i].y**2)
        c=(1+r)/r
        x,y,z = c*straight[i].x,c*straight[i].y,straight[i].z
        straight+=Atom('H',(x,y,z))
#view(straight)
calc = Hotbit(SCC=SCC,txt='chiral.cal',kpts=(1,1,nkpts),gamma_cut=cut)
straight.set_calculator(calc)
e1 = straight.get_potential_energy()



# same thing, but calculate by twisting 2*pi/5 while translating
height = straight.get_cell()[2,2]
chiral = Atoms(container='Chiral')
chiral += straight
chiral.set_container(height=height,angle=2*pi/5)
calc = Hotbit(SCC=SCC,txt='chiral.cal',kpts=(1,1,nkpts),gamma_cut=cut)
chiral.set_calculator(calc)
view(chiral)
e2 = chiral.get_potential_energy()
assert abs(e1-e2)<1E-6


# check the energy conservation for the chiral situation
chiral.rattle(0.1)
calc = Hotbit(SCC=SCC,txt='chiral.cal',width=0.1,kpts=(1,1,1),gamma_cut=cut) #,verbose_SCC=True)
chiral.set_calculator(calc)
assert check_energy_conservation(chiral,dt=0.5*fs,steps=50,tol=0.02,plot=True)
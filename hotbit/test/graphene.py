from ase import *
from hotbit import *
from box.md import check_energy_conservation
from box.systems import graphene
from numpy import *
from hotbit.test.misc import default_param

R=1.416552
nkpts = 10

# graphene's cohesion energy
atoms=graphene(2,2,R)
atoms=Atoms('C2',[(0,0,0),(2*cos(pi/6)*R,R,0)],pbc=(True,True,False),\
            cell=array([[2*cos(pi/6)*R,0,0],[R*cos(pi/6),R*1.5,0],[0,0,5]]) )

calc=Hotbit(SCC=False,txt='graphene.cal',kpts=(nkpts,nkpts,1),**default_param) 
atoms.set_calculator(calc)
coh = atoms.get_potential_energy()/2
assert abs(-9.626283-coh)<1E-6



# energy conservation
atoms[0].z+=0.1
calc=Hotbit(SCC=False,txt='graphene.cal',kpts=(3,3,1),**default_param) 
atoms.set_calculator(calc)
assert check_energy_conservation(atoms,dt=0.5*units.fs,steps=30,tol=1E-2,plot=False)



from ase import *
from hotbit import *
from pylab import *
from numpy import *
from copy import copy

a=4.08
d = a/2
cell=array([[d,d,0],[0,d,d],[d,0,d]])
atoms = Atoms('Au',[(0,0,0)],cell=cell,pbc=(True,True,True))
write('bulk.traj',atoms)

#view(atoms)
calc = Hotbit(SCC=False,kpts=(6,6,6))
atoms.set_calculator(calc)
e=[]
cell0=copy(cell)

for x in linspace(0.9,1.4,10):
    cell = cell0*x
    atoms.set_cell(cell,scale_atoms=True)
    e.append( atoms.get_potential_energy() )
    
plot(linspace(0.9,1.4,10)*a,e)
show()
    

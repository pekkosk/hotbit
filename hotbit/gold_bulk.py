from ase import *
from hotbit import *
from pylab import *

atoms=Atoms('Au',[(0,0,0)],pbc=True)
a = 4.08
b = a/2
atoms.set_cell([[0,b,b],[b,0,b],[b,b,0]])
#view(atoms)

e=[]
kl=[1,2,3,4,6,8,10]
for nk in kl:
    print 'nk',nk
    calc=Hotbit(SCC=False,txt='gold_bulk.cal',kpts=(nk,nk,nk))
    atoms.set_calculator(calc)
    e.append(atoms.get_potential_energy())
    
atoms=Atoms('Au',[(0,0,0)])
calc=Hotbit(SCC=False)
atoms.set_calculator(calc)
e1=atoms.get_potential_energy()

plot(kl,e-e1)
title('Au binding energy in bulk')
xlabel('number of k-points in each direction')
show()


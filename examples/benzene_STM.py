from pylab import *
from ase import *
from hotbit import *
from box.grids import GridData
from ase.data.molecules import molecule

atoms = molecule('C6H6')
calc = Hotbit(SCC=True,width=0.05,txt='benzene.cal')
atoms.set_calculator(calc)
atoms.center(vacuum=3)
atoms.get_potential_energy()
#view(atoms)
r = atoms.get_positions()

calc.set_grid(h=0.4)
ldos = calc.get_grid_LDOS(bias=-6)

current=4 # current in nanoamperes
crit = 2E-4*sqrt(current)
stm = GridData(atoms,data=ldos)  
h = stm.scan(crit,bottom=r[:,2].mean()).transpose()

X,Y,Z = stm.get_grids()
contourf(X,Y,h,50)
hot()
scatter(r[:,0],r[:,1],color='blue',marker='o',s=2)
xlabel(r'x ($\AA$)')
ylabel(r'y ($\AA$)')
xlim(xmin=0,xmax=X[-1])
ylim(ymin=0,ymax=Y[-1])
savefig('benzene_STM.png')


from pylab import *
from ase import *
from hotbit import *

d1,d2 = 1.06, 1.203
atoms = Atoms('HCCH',[(0,0,0),(d1,0,0),(d1+d2,0,0),(2*d1+d2,0,0)])
atoms.center(vacuum=3)
calc = Hotbit(SCC=True,txt='-')
atoms.set_calculator(calc)
atoms.get_potential_energy()

#
# charges and energetics
#
print 'charges:',calc.get_dq()
print 'H: mulliken charge',calc.get_atom_mulliken(0)
print 'C: mulliken charge',calc.get_atom_mulliken(1)
print 'H: atom and bond energy',calc.get_atom_and_bond_energy(0)
print 'C: atom and bond energy',calc.get_atom_and_bond_energy(1)
print 'H: atom energy',calc.get_atom_energy(0)
print 'C: atom energy',calc.get_atom_energy(1)

#
# grid analysis
#
calc.set_grid(h=0.25)
write('pi-bond.cube',atoms,data=calc.get_grid_wf(3))


#
# DOS
#
e,dos,pdos=calc.get_density_of_states(broaden=True,projected=True,width=0.3)
plot(e,dos,label='total DOS')
for l in range(3):
    plot(e,pdos[l,:],label='l=%i' %l)
legend()
savefig('DOS.pdf')
clf()

for l1 in range(2):
    for l2 in range(2):
        x,y = calc.get_covalent_energy(mode='angmom',i=l1,j=l2,width=0.2)
        plot(x,y,label='ecov angmom %i-%i' %(l1,l2))
legend()
savefig('covalent_bonding.pdf')


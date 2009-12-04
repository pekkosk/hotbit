from ase import *
from hotbit import *

atoms = Atoms('H2',[(0,0,0),(1,0,0)])
calc = Hotbit(SCC=True,txt='external.cal')
atoms.set_calculator(calc)

def phi(r,t):
    return 0.78*r[0]

calc.env.add_phi(phi)
#print calc.get_dq(atoms)
calc.solve_ground_state(atoms)
assert calc.get_dq()[1]>0.99

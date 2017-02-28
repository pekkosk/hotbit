from ase import *
from hotbit import Hotbit
from hotbit.test.misc import default_param

calc=Hotbit(SCC=True,txt='test.cal',**default_param)
Au2=Atoms('Au2',positions=[(0,0,0),(1.6,1.2,0.2)],cell=(10,10,10),pbc=False)
Au2.center(vacuum=10)
Au2.set_calculator(calc)
print(Au2.get_potential_energy())

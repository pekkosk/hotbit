from ase import *
from box import Atoms
from box import mix
from hotbit import Calculator
#from hotbit import Calculator0
from ase.units import Bohr, Hartree
from hotbit.test.misc import default_param

# C-non-SCC
calc=Calculator(SCC=False,txt='test.cal',**default_param)
C=Atoms('C',positions=[(0,0,0)],cell=(10,10,10),pbc=False)
C.center(vacuum=100)
C.set_calculator(calc)
e=C.get_potential_energy()
if abs(-37.8891472348-e)>1E-6:
    raise RuntimeError('energy %f, should be %f'  %(e,-37.8891472348))
    
# C, SCC
calc=Calculator(SCC=True,txt='test.cal',**default_param)
C=Atoms('C',positions=[(0,0,0)],cell=(10,10,10),pbc=False)
C.center(vacuum=100)
C.set_calculator(calc)
e=C.get_potential_energy()
if abs(-37.8891472348-e)>1E-6:
    raise RuntimeError('energy %f, should be %f'  %(e,-37.8891472348))


# rotate Au-dimer
#calc0=Calculator0(SCC=True,txt='test.cal',**default_param)
calc=Calculator(SCC=True,txt='test.cal',**default_param)
Au2=Atoms('Au2',positions=[(0,0,0),(2.6,0,0)],cell=(10,10,10),pbc=False)
Au2.center(vacuum=10)
#Au2.set_calculator(calc0)
e=-152.981553763
for i in range(10):    
    Au2.set_calculator(calc)
    vector=mix.random_direction(vector=True)
    Au2.rotate(vector*random.rand()*2*pi)
    Au2.center(vacuum=10)
    e2=Au2.get_potential_energy()
    if abs(e-e2)>1E-4:
        raise RuntimeError('energy for Au2 %f, should be %f (while rotating)' %(e2,e))
    

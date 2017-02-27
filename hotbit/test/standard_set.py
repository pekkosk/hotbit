from hotbit import Hotbit 
#from hotbit import Calculator0
from ase import *
from ase.structure import molecule
from hotbit.test.misc import default_param
import sys

systems=['C2H6','CH3CH2O','H2COH','H2','C2H6CHOH','isobutane']
energies0=[-153.471502365,-232.206844117,-166.801950557,-19.356876265,-313.487801685,-287.534640343]
free = [-113.81982575044,-191.66618644170001,-141.09652875422,-12.680510375480001,-248.57609931692002,-214.95914112540004]
energies = [e-ef for e,ef in zip(energies0,free)]

eps=1E-3
for system,e in zip(systems,energies):
    atoms=molecule(system)
    atoms.center(vacuum=10)
    atoms.set_pbc(False)    
    if e==0:
        calc=Calculator0(verbose=True,SCC=False,txt='standard.cal',**default_param)
        atoms.set_calculator(calc)
        print('new system',system,atoms.get_potential_energy())
        sys.exit(0)        
        
    calc=Hotbit(verbose=True,SCC=True,txt='standard.cal',**default_param)
    atoms.set_calculator(calc)
    e1=atoms.get_potential_energy()
    if abs(e1-e)>eps:
        raise AssertionError('Energy for %s is %.7f, while should be %.7f' %(system,e1,e))


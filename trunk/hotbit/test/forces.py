from hotbit import Calculator
#from hotbit import Calculator0
from ase import *
from hotbit.test.misc import default_param
from hotbit.test.misc import molecule
from numpy import *
from box.md import microcanonical, energy_conservation
from ase.units import fs
import sys
           
systems=['H2COH','AuC']          

#default_param['convergence']=1E-5
#default_param['Anderson_memory']=3
#default_param['width']=0.01

for SCC in [False,True]:            
    for system in systems:   
        print '    ... forces for %s, SCC=' %system, SCC         
        calc=Calculator(verbose=True,SCC=True,txt='forces.cal',**default_param)
        atoms=molecule(system)
        atoms.center(vacuum=5)
        atoms[0].z+=0.2
        atoms=Atoms(atoms)
        atoms.set_calculator(calc)
        sys.stdout.flush()
        rec, de, du=energy_conservation(atoms,dt=0.2*fs,steps=50)
        calc.__del__()
        assert de/du<0.01
        
    





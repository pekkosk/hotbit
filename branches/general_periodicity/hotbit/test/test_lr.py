from hotbit import Calculator
from hotbit import Calculator0
from hotbit import LinearResponse
from hotbit.test.misc import default_param
from ase import *
from box import Atoms

import numpy as nu
from box import mix
import pylab as pl
from ase.units import Hartree
 
Na3=Atoms('Na3',positions=[(1.69649997,0,0),(-1.69649997,0,0),(0,2.9384241,0)],cell=(50,50,50),pbc=False)  
tm=mix.Timer()

calc=Calculator(charge=1,SCC=True,txt='test_lr.cal',**default_param)
Na3.set_calculator(calc)
calc.solve_ground_state(Na3)
lr=LinearResponse(calc)
omega,F=lr.get_linear_response()
e,f=mix.broaden(omega,F,width=0.1)

pl.scatter(e2,f2)
pl.plot(e,f,label='python')
pl.legend()
pl.show()
       
calc=Calculator(SCC=True,txt='test_lr.cal',**default_param)
C60=Atoms(read('C60.xyz'))
C60.set_calculator(calc)
calc.solve_ground_state(C60)
lr=LinearResponse(calc,energy_cut=10/Hartree)
omega,F=lr.get_linear_response()
e,f=mix.broaden(omega,F,width=0.1)
pl.plot(e,f)
pl.show()


 




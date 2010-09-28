from hotbit import Hotbit
from hotbit import LinearResponse
from hotbit.test.misc import default_param
from ase import *
from box import Atoms
plot=True
try:
    import pylab as pl
except:
    plot=False
    
import numpy as np
from box import mix

    
from ase.units import Hartree
 
Na3=Atoms('Na3',positions=[(1.69649997,0,0),(-1.69649997,0,0),(0,2.9384241,0)],cell=(50,50,50),pbc=False)  
tm=mix.Timer()

calc=Hotbit(charge=1,SCC=True,txt='test_lr.cal',**default_param)
Na3.set_calculator(calc)
calc.solve_ground_state(Na3)
lr=LinearResponse(calc)
omega,F=lr.get_linear_response()
e,f=mix.broaden(omega,F,width=0.1)

if plot:
    pl.scatter(e2,f2)
    pl.plot(e,f,label='python')
    pl.legend()
    pl.show()
       
calc=Hotbit(SCC=True,txt='test_lr.cal',**default_param)
C60=Atoms(read('C60.xyz'))
C60.set_calculator(calc)
calc.solve_ground_state(C60)
lr=LinearResponse(calc,energy_cut=10/Hartree)
omega,F=lr.get_linear_response()
e,f=mix.broaden(omega,F,width=0.1)

if plot:
    pl.plot(e,f)
    pl.show()


 




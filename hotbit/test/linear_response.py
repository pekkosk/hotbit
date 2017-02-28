from ase import *
from hotbit import Hotbit
from hotbit.analysis import LinearResponse
from hotbit.test.misc import default_param

atoms=Atoms('Na3',[(1.6964999745231999,0,0),(-1.6964999745231999,0,0),(0,2.9384240999630005,0)])
atoms.center(vacuum=5)

default_param['width'] = 0.0136

calc=Hotbit(SCC=True,charge=1,txt='linear_response.cal',**default_param)
atoms.set_calculator()
calc.solve_ground_state(atoms)

lr=LinearResponse(calc,energy_cut=2000,txt='linear_response.txt')
lr.run()
lr.plot_spectrum('Na3+_lr.png',width=0.08)

el = [1.81951,1.81951,7.43599,7.43599]
fl = [0.43036,0.43036,4.27744,4.27744]

for i in range(4):
    e,f = lr.get_excitation(i)
    assert abs(e-el[i])<1E-4 and abs(f-fl[i])<1E-4
    

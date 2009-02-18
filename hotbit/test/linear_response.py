from ase import *
from hotbit import Hotbit
from hotbit.analysis import LinearResponse
from ase.data.molecules import molecule
from hotbit.test.misc import default_param

atoms=Atoms('Na3',[(1.6964999745231999,0,0),(-1.6964999745231999,0,0),(0,2.9384240999630005,0)])
atoms.center(vacuum=5)

default_param['width'] = 0.0136

calc=Hotbit(SCC=True,charge=1,**default_param)
atoms.set_calculator()
calc.solve_ground_state(atoms)

lr=LinearResponse(calc,energy_cut=2000)
lr.run()
lr.plot_spectrum('Na3+_lr.png',width=0.08)
lr.info()

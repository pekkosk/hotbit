"""
Test the DFTB file format reader.
"""

import numpy as np

import ase
from ase.units import Hartree
from hotbit import Hotbit

###


dir1 = '/Users/pas/Sourcen/hotbit/param/'
#dir2 = '/Users/pas/Sourcen/hotbit/param/'
dir2 = '/Users/pas/Sourcen/mdcore/data/slater_koster_tables/Frauenheim/download-2009-10-23/mio-0-1/'

###

a = ase.molecule('C6H6')
a.center(vacuum=6.0)

c = Hotbit(
    elements = {
        'H': dir2 + 'H-H.skf',
        'C': dir2 + 'C-C.skf'
        },
    tables   = {
        'HH': dir2 + 'H-H.skf',
        'CH': dir2 + 'C-H.skf',
        'CC': dir2 + 'C-C.skf'
        },
    SCC      = False
    )
a.set_calculator(c)

ase.FIRE(a).run(fmax=0.01)

ase.write('C6H6.cfg', a)

###

for i in c.el.get_present():
    el = c.el.get_element(i)
    sym = el.get_symbol()
    orbs = el.get_valence_orbitals()
    for orb in orbs:
        print(sym, orb, el.get_epsilon(orb), el.get_epsilon(orb)*Hartree)

x = np.linspace(1.0, 12.0, 1000)

y, dy = c.ia.h['CC'](x)
d = np.transpose(np.append(x.reshape(1, -1), y, axis=0))
np.savetxt('CC_H.out', d)

y, dy = c.ia.s['CC'](x)
d = np.transpose(np.append(x.reshape(1, -1), y, axis=0))
np.savetxt('CC_S.out', d)

y = c.rep.vrep['CC'](x)
d = np.transpose([x, y])
np.savetxt('CC_rep.out', d)

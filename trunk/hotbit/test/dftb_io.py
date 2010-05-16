"""
Test the DFTB file format reader.
"""

import ase
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
        'H': dir1 + 'H.elm',
        'C': dir1 + 'C.elm'
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

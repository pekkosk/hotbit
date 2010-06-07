"""
Test the mio parametrization of Frauenheim and co-workers.
"""

import sys
import glob

import numpy as np

from ase import read, FIRE, QuasiNewton
from hotbit import Hotbit

###

# From Elstner et al., Phys. Rev. B 58, 7260
noscc_db = {
    'C=O': 1.296,
    'C-N': 1.296,
    'N-H': 1.003,
    'C-H': 1.130,
    'OCN': 127.0
    }

scc_db = {
    'C=O': 1.224,
    'C-N': 1.382,
    'N-H': 0.996,
    'C-H': 1.131,
    'OCN': 125.5
    }

db = {
    False: noscc_db,
    True: scc_db
    }

def check_q(db, name, value):
    refvalue = db[name]
    print '%10s  %10.3f  %10.3f' % ( name, value, refvalue )
    #assert abs(value-refvalue) < 1e-3

###

def params_for_database_from(path):
    fns = glob.glob("%s/*-*.skf" % path)

    elements = { }
    tables = { }
    for fn in fns:
        i0 = fn.rfind('/')
        i1 = fn.rfind('-')
        i2 = fn.rfind('.')
        el1 = fn[i0+1:i1]
        el2 = fn[i1+1:i2]

        if el1 == el2:
            elements[el1] = fn
        tables['%s%s' % ( el1, el2 )] = fn

    return { 'elements': elements, 'tables': tables }

###

if len(sys.argv) != 2:
    print "Syntax: mio.py <path to database>"
    sys.exit(999)

params = params_for_database_from(sys.argv[1])


for SCC in [ False, True ]:
    if SCC:
        print "--- SCC ---"
    else:
        print "--- no SCC ---"

    calc = Hotbit(
        charge_density = 'Slater',
        SCC = SCC,
        #verbose = True,
        txt = 'mio.out',
        **params)

    a = read('formamide.xyz')
    a.center(vacuum=10.0)
    a.set_pbc(False)
    a.set_calculator(calc)

    QuasiNewton(a, logfile='fire.log').run(fmax=0.01)

    iO = 0 
    iC = 1
    iN = 2
    iHC = 3
    iHN = 4

    assert a[iO].get_symbol() == 'O'
    assert a[iC].get_symbol() == 'C'
    assert a[iN].get_symbol() == 'N'
    assert a[iHC].get_symbol() == 'H'
    assert a[iHN].get_symbol() == 'H'

    check_q(db[SCC], 'C=O', a.get_distance(iC, iO))
    check_q(db[SCC], 'C-N', a.get_distance(iC, iN))
    check_q(db[SCC], 'N-H', a.get_distance(iN, iHN))
    check_q(db[SCC], 'C-H', a.get_distance(iC, iHC))


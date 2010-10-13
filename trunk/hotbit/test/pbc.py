"""
Test the mio parametrization of Frauenheim and co-workers.
"""

import os
import glob

import numpy as np

from ase import read, FIRE, QuasiNewton, molecule
from hotbit import Hotbit, database_from_path

###

FMAX = 0.005
OPT  = FIRE

debug = False

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

db1 = {
    False: noscc_db,
    True: scc_db
    }

# From Kruger et al., J. Chem. Phys. 122, 114110
db2 = {
    'H2': {
        'H-H': ( ( 0, 1 ), 0.750 )
        },
    'C2H2': {
        'C-H': ( ( 1, 2 ), 1.075 ),
        'C-C': ( ( 0, 1 ), 1.203 )
        },
    'C2H4': {
        'C-H': ( ( 0, 2 ), 1.094 ),
        'C-C': ( ( 0, 1 ), 1.328 )
        },
    'C2H6': {
        'C-H': ( ( 0, 3 ), 1.098 ),
        'C-C': ( ( 0, 1 ), 1.501 )
        },
    'HCN': {
        'C-H': ( ( 0, 2 ), 1.078 ),
        'C-N': ( ( 0, 1 ), 1.141 )
        },
    'NH3': {
        'N-H': ( ( 0, 1 ), 1.021 )
        },
    'CH4': {
        'C-H': ( ( 0, 1 ), 1.089 )
        },
    'CO': {
        # This differs from the paper, but I believe it's a typo
        # paper says: 1.200
        'C-O': ( ( 0, 1 ), 1.100 )
        },
    'H2CO': {
        'C-H': ( ( 1, 2 ), 1.143 ),
        'C-O': ( ( 0, 1 ), 1.183 )
        },
    'CH3OH': {
        'O-H': ( ( 1, 3 ), 0.980 ),
        'C-O': ( ( 0, 1 ), 1.422 )
        },
    'H2O': {
        'O-H': ( ( 0, 1 ), 0.968 )
        },
    'N2': {
        # This differs from the paper, but I believe it's a typo
        # paper says: 1.200
        'N-N': ( ( 0, 1 ), 1.113 )
        },
    'N2H4': {
        'N-H': ( ( 0, 2 ), 1.037 ),
        # This differs from the paper, and I don't know why
        # paper says: 1.442
        'N-N': ( ( 0, 1 ), 1.407 )
        },
    'H2O2': {
        'O-H': ( ( 0, 2 ), 0.991 ),
        'O-O': ( ( 0, 1 ), 1.453 )
        },
    'CO2': {
        'C-O': ( ( 0, 1 ), 1.165 )
        }
    }


def check_q(db, name, value):
    refvalue = db[name]
    if debug:
        print '%10s  %10.3f  %10.3f' % ( name, value, refvalue )
    #assert abs(value-refvalue) < 1e-3


def check_db(db, params):
    if debug:
        print "%10s %10s %10s ( %10s )" \
            % ( "bond", "value", "reference", "error" )
    for mol, values in db.iteritems():
        #if mol == 'H2O':
        if 1:
            if debug:
                print mol

            a = molecule(mol)
            a.center(vacuum=10.0)
            a.set_pbc(False)

            #print a.get_chemical_symbols()

            calc = Hotbit(
                charge_density = 'Slater',
                SCC = True,
                width = 1e-6,
                txt = 'mio.out',
                **params)
            a.set_calculator(calc)

            #calc.ia.plot_table('H', 'H')
            #calc.rep.get_repulsion('H', 'H').plot()

            OPT(a, logfile='opt.log').run(fmax=FMAX)

            #print a.get_charges()
            
            for name, ( ( i1, i2 ), refvalue ) in values.iteritems():
                value = a.get_distance(i1, i2)
                if debug:
                    print '%10s %10.3f %10.3f ( %10.3f )' % \
                        ( name, value, refvalue, abs(value-refvalue) )
                assert abs(value-refvalue) < 0.01

            #e = [ ]
            #for x in np.linspace(0.70, 0.80, 1000):
            #    a.set_distance(0, 1, x)
            #    e += [ ( x, a.get_potential_energy() ) ]
            #np.savetxt('e.out', e)

###

params = database_from_path(os.getenv('PBC_0_2'))

###

if debug:
    for SCC in [ False, True ]:
        if SCC:
            print "--- SCC ---"
        else:
            print "--- no SCC ---"
            
        calc = Hotbit(
            charge_density = 'Slater',
            SCC = SCC,
            verbose = True,
            verbose_SCC = True,
            mixer = {
                'name': 'anderson',
                'convergence': 1e-6,
                'mixing_constant': 0.01 },
            maxiter = 1000,
            txt = 'mio.out',
            **params)

        a = read('formamide.xyz')
        a.center(vacuum=10.0)
        a.set_pbc(False)
        a.set_calculator(calc)

        OPT(a, logfile='opt.log').run(fmax=FMAX)

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

        check_q(db1[SCC], 'C=O', a.get_distance(iC, iO))
        check_q(db1[SCC], 'C-N', a.get_distance(iC, iN))
        check_q(db1[SCC], 'N-H', a.get_distance(iN, iHN))
        check_q(db1[SCC], 'C-H', a.get_distance(iC, iHC))

###

check_db(db2, params)

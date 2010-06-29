#! /usr/bin/env python

from math import sqrt

import numpy as np

from ase import Atoms, write
from ase.lattice.compounds import CsCl, NaCl, ZnS

from hotbit.atoms import Atoms as HotbitAtoms
from hotbit.coulomb import EwaldSum, MultipoleExpansion

###

L_MAX  = 8
K      = 5

Q      = 1.0
a0     = 1.0

debug  = False

###

# FIXME!!! Look for more digits
M_NaCl  = 1.747565
M_CsCl  = 1.762675
M_ZnS   = 1.638055

systems = [
#    ( "NaCl", M_NaCl, 0.5,
#       NaCl(['Na', 'Cl'],
#            latticeconstant = a0,
#            size            = [1, 1, 1]) ),
    ( "CsCl", M_CsCl, sqrt(3.0)/2,
      CsCl(['Cs', 'Cl'],
          latticeconstant   = a0,
          size              = [1, 1, 1]) ),
#    ( "ZnS", M_ZnS, sqrt(3.0)/4,
#      ZnS(['Zn', 'S'],
#          latticeconstant   = a0,
#          size              = [1, 1, 1]) ),
# Some tests with extended unit cells, not necessary
#    ( "large NaCl", M_NaCl, 0.5,
#       NaCl(['Na', 'Cl'],
#            latticeconstant = a0,
#            size            = [4, 4, 4]) ),
#    ( "large CsCl", M_CsCl, sqrt(3.0)/2,
#      CsCl(['Cs', 'Cl'],
#          latticeconstant   = a0,
#          size              = [4, 4, 4]) ),
#    ( "large ZnS", M_ZnS, sqrt(3.0)/4,
#      ZnS(['Zn', 'S'],
#          latticeconstant   = a0,
#          size              = [4, 4, 4]) )
      ]

solvers = [
    MultipoleExpansion(L_MAX, 3, K),
#    MultipoleExpansion(L_MAX, 4, K),
    EwaldSum(12, 0.001)
]

if debug:
    print "%20s   %8s  %8s  (%8s)" % ( "compound", "M", "ref.", "error" )
    print "%20s   %8s  %8s  %10s" % ( "-------------------",
                                       "--------", "--------", "----------" )

for sol in solvers:
    if debug:
        print "=== %s ===" % sol.__class__
    for name, target_M, nnd, a in systems:
        syms = a.get_chemical_symbols()
        
        #a.set_charges([ (Q if sym == syms[0] else -Q) for sym in syms ])
        # to work with older versions
        a.set_charges([ (-Q,Q)[sym==syms[0]] for sym in syms ])

        a.translate([0.25*a0,0.25*a0,0.25*a0])
        if debug:
            write("%s.cfg" % name, a)

        ha = HotbitAtoms(a, container='Bravais')

        sol.update(ha, ha.get_charges())

        phi  = sol.get_potential()
        e    = np.sum(a.get_charges()*phi)/2
        M    = -2*e*a0*nnd/(len(a))
        err  = abs(M-target_M)

        if debug:
            print "%20s   %8.6f  %8.6f  (%8.6e)" % ( name, M, target_M, err )

        assert err < 1e-3

            

if debug:
    for sol in solvers:
        sol.timer.summary()

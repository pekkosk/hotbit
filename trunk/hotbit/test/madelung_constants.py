#! /usr/bin/env python

from math import sqrt

import numpy as np

from ase import Atoms, write
from ase.units import Hartree, Bohr
from ase.lattice.compounds import CsCl, NaCl, ZnS

from hotbit.atoms import Atoms as HotbitAtoms
from hotbit.multipole_expansion import MultipoleExpansion

from hotbit.ewald_sum import EwaldSum


###

L_MAX  = 8
N      = 5
K      = 5

Q      = 1.0
a0     = 1.0

debug  = True

###

# FIXME!!! Look for more digits
M_NaCl  = 1.747565
M_CsCl  = 1.762675
M_ZnS   = 1.6381

systems = [
    ( "NaCl", M_NaCl, 0.5,
       NaCl(['Na', 'Cl'],
            latticeconstant = a0,
            size            = [1, 1, 1]) ),
    ( "CsCl", M_CsCl, sqrt(3.0)/2,
      CsCl(['Cs', 'Cl'],
          latticeconstant   = a0,
          size              = [1, 1, 1]) ),
    ( "ZnS", M_ZnS, sqrt(3.0)/4,
      ZnS(['Zn', 'S'],
          latticeconstant   = a0,
          size              = [1, 1, 1]) ),
    ( "large NaCl", M_NaCl, 0.5,
       NaCl(['Na', 'Cl'],
            latticeconstant = a0,
            size            = [4, 4, 4]) ),
    ( "large CsCl", M_CsCl, sqrt(3.0)/2,
      CsCl(['Cs', 'Cl'],
          latticeconstant   = a0,
          size              = [4, 4, 4]) ),
    ( "large ZnS", M_ZnS, sqrt(3.0)/4,
      ZnS(['Zn', 'S'],
          latticeconstant   = a0,
          size              = [4, 4, 4]) )
      ]

solvers = [
    MultipoleExpansion(L_MAX, N, K),
    EwaldSum(12, 0.001)
]

if debug:
    print "%20s   %8s  %8s  (%8s)" % ( "compound", "M", "ref.", "error" )
    print "%20s   %8s  %8s  %10s" % ( "-------------------",
                                       "--------", "--------", "----------" )

for sol in solvers:
    print "=== %s ===" % sol.__class__
    for name, target_M, nnd, a in systems:
        syms = a.get_chemical_symbols()
        a.set_charges([ Q if sym == syms[0] else -Q for sym in syms ])

        a.translate([0.25*a0,0.25*a0,0.25*a0])
        if debug:
            write("%s.cfg" % name, a)

        ha = HotbitAtoms(a, container='Bravais')

        sol.update(ha, ha.get_charges())

        phi  = sol.get_potential()
        e    = np.sum(a.get_charges()*phi)
        M    = -e*a0*nnd/(len(a)*Hartree*Bohr)
        err  = abs(M-target_M)

        if debug:
            print "%20s   %8.6f  %8.6f  (%8.6e)" % ( name, M, target_M, err )

        assert err < 1e-3

            

if debug:
    for sol in solvers:
        sol.timer.summary()

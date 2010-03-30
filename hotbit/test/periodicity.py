#! /usr/bin/env python

import numpy as np

from ase import Atoms
from hotbit.atoms import Atoms as HotbitAtoms
from hotbit.multipole import get_moments
from hotbit.multipole_periodicity import MultipolePeriodicity

###

NAT         = 3
SX, SY, SZ  = 10.0, 10.0, 10.0
CHARGE      = 1.0

L_MAX       = 8
N           = 5
K           = 5

THRES_MOM   = 1e-6
THRES_PHI   = 1e-4
THRES_E     = 1e-4

###

def bravais_test(a):
    b = HotbitAtoms(
        atoms      = a,
        container  = 'Bravais'
        )

    # Check multipole moments
    mp  = MultipolePeriodicity(L_MAX, 3, 3)
    mp.update(b, b.get_charges())

    moments          = mp.get_moments()
    M0_l_mp, M_L_mp  = moments[1]

    # Next neighbor shell by multipole expansion
    mp  = MultipolePeriodicity(L_MAX, 3, 2)
    mp.update(b, b.get_charges())

    phi_mp, E_mp     = mp.get_potential_and_field()

    # Next neighbor shell by direct summation
    mp  = MultipolePeriodicity(L_MAX, 9, 1)
    mp.update(b, b.get_charges())

    phi_dir, E_dir  = mp.get_potential_and_field()

    # Multipole moments from large cell
    rep = [ 3 if i else 1 for i in a.get_pbc() ]
    a *= rep

    M0_l, M_L  = get_moments(a.get_positions(), a.get_charges(), L_MAX, np.sum(a.get_cell(), axis=0)/2)

    assert np.max(np.abs(M0_l-M0_l_mp)) < THRES_MOM
    assert np.max(np.abs(M_L-M_L_mp)) < THRES_MOM

    # Compare fields and potentials obtained by the multipole expansion
    # and from direct summation

    err_phi = np.max(np.abs(phi_mp-phi_dir))
    err_E   = np.max(np.abs(E_mp-E_dir))

    print "error(phi)  = ", err_phi
    print "error(E)    = ", err_E

    assert err_phi < THRES_PHI
    assert err_E < THRES_E



q  = (2*np.random.random([NAT])-1)*CHARGE
q -= np.sum(q)/len(q)

# 1D periodicity
a = Atoms('%iH' % NAT,
          positions  = np.random.random([NAT,3])*SX,
          charges    = q,
          cell       = [ SX, SY, SZ ],
          pbc        = [ False, False, True ]
          )

bravais_test(a)


# 2D periodicity
a = Atoms('%iH' % NAT,
          positions  = np.random.random([NAT,3])*SX,
          charges    = q,
          cell       = [ SX, SY, SZ ],
          pbc        = [ True, False, True ]
          )

bravais_test(a)


# 3D periodicity
a = Atoms('%iH' % NAT,
          positions  = np.random.random([NAT,3])*SX,
          charges    = q,
          cell       = [ SX, SY, SZ ],
          pbc        = True
          )

bravais_test(a)



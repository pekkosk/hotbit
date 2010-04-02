#! /usr/bin/env python

from math import pi

import numpy as np

from ase import Atoms, write
from hotbit.atoms import Atoms as HotbitAtoms
from hotbit.multipole import get_moments
from hotbit.multipole_expansion import MultipoleExpansion

from box.fd_forces import check_forces

###

NAT         = 4
SX, SY, SZ  = 10.0, 10.0, 10.0
CHARGE      = 1.0

L_MAX       = 8

NRUNS       = 1

debug  = False

###

TOL_MOM = {
    3: 1e-5,
    5: 1e-7
    }
TOL_PHI = {
    3: 1e-4,
    5: 1e-6
    }
TOL_E = {
    3: 1e-4,
    5: 1e-6
    }
TOL_FOR = 1e-6

###

def electrostatics_test(b, r=3, r0=None):
    if r0 is None:
        r0 = np.sum(b.get_cell(), axis=0)/2

    # Check multipole moments
    mp  = MultipoleExpansion(L_MAX, r, 3)
    mp.update(b, b.get_charges(), r0)

    # Store moments for later
    moments          = mp.get_moments()
    M0_l_mp, M_L_mp  = moments[1]

    # Check if the field is the derivative of the potential
    b.set_calculator(mp)
    ffd, f0, err = check_forces(b)
    if debug:
        print "Finite differences forces:"
        print ffd
        print "Analytical forces:"
        print f0
        print "Error:"
        print err

    assert err < TOL_FOR

    # Next neighbor shell by multipole expansion
    mp  = MultipoleExpansion(L_MAX, r, 2)
    mp.update(b, b.get_charges(), r0)

    phi_mp, E_mp     = mp.get_potential_and_field()

    # Next neighbor shell by direct summation
    mp  = MultipoleExpansion(L_MAX, 3*r, 1)
    mp.update(b, b.get_charges(), r0)

    phi_dir, E_dir  = mp.get_potential_and_field()

    # Multipole moments from large cell,
    # transform symmetrically around the origin
    rep  = [ (r-1)/2 if i else 0 for i in b.get_pbc() ]
    rep  = [ (-(r-1)/2, (r-1)/2) if i else (0, 0) for i in b.get_pbc() ]
    c    = b.extended_copy(tuple(rep))

    M0_l, M_L  = get_moments(c.get_positions(), c.get_charges(), L_MAX, r0)

    #print M_L
    #print M_L_mp

    err_mom0  = np.max(np.abs(M0_l-M0_l_mp))
    err_mom   = np.max(np.abs(M_L-M_L_mp))

    if debug:
        print "error(mom)  = ", err_mom0, err_mom

    assert err_mom0 < TOL_MOM[r]
    assert err_mom < TOL_MOM[r]

    # Compare fields and potentials obtained by the multipole expansion
    # and from direct summation

    err_phi  = np.max(np.abs(phi_mp-phi_dir))
    err_E    = np.max(np.abs(E_mp-E_dir))

    if debug:
        print "error(phi)  = ", err_phi
        print "error(E)    = ", err_E

    assert err_phi < TOL_PHI[r]
    assert err_E < TOL_E[r]



for i in range(NRUNS):
    q  = (2*np.random.random([NAT])-1)*CHARGE
    q -= np.sum(q)/len(q)

    if True:
        # no periodicity, this should simply not fail
        if debug:
            print "0D"
        a = Atoms('%iH' % NAT,
                  positions  = np.random.random([NAT,3])*SX,
                  charges    = q,
                  cell       = [ SX, SY, SZ ],
                  pbc        = False
                  )
        b = HotbitAtoms(
            atoms      = a,
            container  = 'Bravais'
            )

        electrostatics_test(b, 3)
        electrostatics_test(b, 5)


    if True:
        # 1D periodicity
        if debug:
            print "1D"
        a = Atoms('%iH' % NAT,
                  positions  = np.random.random([NAT,3])*SX,
                  charges    = q,
                  cell       = [ SX, SY, SZ ],
                  pbc        = [ False, False, True ]
                  )
        b = HotbitAtoms(
            atoms      = a,
            container  = 'Bravais'
            )

        electrostatics_test(b, 3)
        electrostatics_test(b, 5)


    if True:
        # 1D and twisted periodicity
        if debug:
            print "1D - twisted"
        a = Atoms('%iH' % NAT,
                  positions  = (np.random.random([NAT,3])-0.5)*SX,
                  charges    = q,
                  cell       = [ SX, SY, SZ ],
                  pbc        = [ False, False, True ]
                  )
        b = HotbitAtoms(
            atoms      = a,
            container  = 'Chiral'
            )
        b.set_container(
            angle      = 2*pi/30,
            height     = SZ
            )

        electrostatics_test(b, 3, r0=np.zeros(3))
        electrostatics_test(b, 5, r0=np.zeros(3))


    if True:
        # 2D periodicity
        if debug:
            print "2D"
        a = Atoms('%iH' % NAT,
                  positions  = np.random.random([NAT,3])*SX,
                  charges    = q,
                  cell       = [ SX, SY, SZ ],
                  pbc        = [ True, False, True ]
                  )
        b = HotbitAtoms(
            atoms      = a,
            container  = 'Bravais'
            )

        electrostatics_test(b, 3)
        electrostatics_test(b, 5)


    if False:
        # 3D periodicity
        if debug:
            print "3D"
        a = Atoms('%iH' % NAT,
                  positions  = np.random.random([NAT,3])*SX,
                  charges    = q,
                  cell       = [ SX, SY, SZ ],
                  pbc        = True
                  )
        b = HotbitAtoms(
            atoms      = a,
            container  = 'Bravais'
            )

        electrostatics_test(a, 3)
        electrostatics_test(a, 5)



#! /usr/bin/env python

from math import pi

import numpy as np

from ase import Atoms
from ase.io import  write
from hotbit.atoms import Atoms as HotbitAtoms
from hotbit.coulomb import MultipoleExpansion
from hotbit.coulomb.multipole import get_moments

from box.fd_forces import check_forces, check_field

###

NAT         = 4
SX, SY, SZ  = 10.0, 10.0, 10.0
CHARGE      = 1.0

L_MAX       = 8

NRUNS       = 1

debug  = False

###

TOL_MOM = {
    3: 1e-7,
    5: 1e-7
    }
TOL_PHI = {
    3: 1e-5,
    5: 1e-6
    }
TOL_E = {
    3: 1e-5,
    5: 1e-7
    }
TOL_FOR = 1e-6
#TOL_FIELD = 1e-3

###

def electrostatics_test(b, r=3, r0=None):
    if r0 is None:
        r0 = np.sum(b.get_positions(), axis=0)/len(a)

    k1 = [ (1,1)[i] for i in b.get_pbc() ]
    k3 = [ (1,3)[i] for i in b.get_pbc() ] 

    # Check multipole moments
    mp  = MultipoleExpansion(L_MAX, r, k3, r0)
    mp.update(b, b.get_initial_charges())

    # Store moments and field for later
    moments          = mp.get_moments()
    M0_l_mp, M_L_mp  = moments[1]

    phi_mp, E_mp     = mp.get_potential_and_field(b)

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

#    ffd, f0, err = check_field(b, mp)
#    if debug:
#        print "Finite differences field:"
#        print ffd
#        print "Analytical field:"
#        print f0
#        print "Error:"
#        print err
#
#    assert err < TOL_FIELD

    # Next neighbor shell by direct summation
    mp  = MultipoleExpansion(L_MAX, r*r*r, k1, r0)
    mp.update(b, b.get_initial_charges())

    phi_dir, E_dir  = mp.get_potential_and_field(b)

    # Multipole moments from large cell,
    # transform symmetrically around the origin
    #rep  = [ (r-1)/2 if i else 0 for i in b.get_pbc() ]
    #rep  = [ (-(r-1)/2, (r-1)/2) if i else (0, 0) for i in b.get_pbc() ]
    rep  = [ (0,(r-1)/2)[i] for i in b.get_pbc() ]
    rep  = [ ((0,0),(-(r-1)/2, (r-1)/2))[i] for i in b.get_pbc() ]
    c    = b.extended_copy(tuple(rep))

    M0_l, M_L  = get_moments(c.get_positions(), c.get_initial_charges(), L_MAX, r0)

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
#    r  = np.random.random([NAT,3])*SX
#
#    q  = (2*np.random.random([NAT])-1)*CHARGE
#    q -= np.sum(q)/len(q)

    r  = [ [   SX/4,   SY/4,   SZ/4 ],
           [ 3*SX/4,   SY/4,   SZ/4 ],
           [   SX/4, 3*SY/4,   SZ/4 ],
           [   SX/4,   SY/4, 3*SZ/4 ] ]
    q  = [ 1, -1, 1, -1 ]

    if False:
        # no periodicity, this should simply not fail
        if debug:
            print "0D"
        a = Atoms('%iH' % NAT,
                  positions  = r,
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
                  positions  = r,
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
                  positions  = r,
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

        electrostatics_test(b, 3)
        electrostatics_test(b, 5)


    if True:
        # 2D periodicity
        if debug:
            print "2D"
        a = Atoms('%iH' % NAT,
                  positions  = r,
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
                  positions  = r,
                  charges    = q,
                  cell       = [ SX, SY, SZ ],
                  pbc        = True
                  )
        b = HotbitAtoms(
            atoms      = a,
            container  = 'Bravais'
            )

        electrostatics_test(b, 3)
        electrostatics_test(b, 5)



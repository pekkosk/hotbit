#
# Test the multipole transformation rules
#

import sys

from math import exp, sqrt, sin, cos, pi, factorial
from cmath import exp as cexp

import random

import numpy as np
import ase

from hotbit.multipole import zero_moments, get_moments

from _hotbit import multipole_to_multipole, multipole_to_local, local_to_local, transform_multipole

###

def D(l, m, n, a, b, g):
    """
    This is Wigner's formula. For debugging purposes only.
    """

    h = 0.0
    for s in range(0, max(l+abs(m), l+abs(n))+1):
        f1 = l+n-s
        f2 = m-n+s
        f3 = l-m-s
        if f1 >= 0 and f2 >= 0 and f3 >= 0:
            h += (-1.0)**(m-n-s)/(factorial(f1)*factorial(s)*factorial(f2)*factorial(f3)) * (cos(b/2))**(2*l+n-m-2*s) * (sin(b/2))**(m-n+2*s)

    return cexp(-1j*m*a) * cexp(-1j*n*g) * h * factorial(l+n)*factorial(l-n)

###

ANGLE       = 10*pi/180

NRUNS       = 10
NAT         = 20
SX, SY, SZ  = 10.0, 10.0, 10.0
CHARGE      = 1.0
L_MAX       = 5

THRES_MOM   = 1e-6
#THRES_PHI   = 1e-3
#THRES_PHI2  = 1e-3
THRES_PHI   = 1
THRES_PHI2  = 1
THRES_ROT   = 1e-9

###

for run in range(NRUNS):
    r0  = np.array( [ SX/2, SY/2, SZ/2 ] )
    r0c = np.array( [ SX, SY, SZ ] )

    # Random atoms and charges (charges between -1 and 1)
    a = [ ]
    for i in range(8):
        a += [ ase.Atoms(
                "%iH" % NAT,
                positions  = np.random.random([NAT,3])*SX,
                charges    = (2*np.random.random([NAT])-1)*CHARGE,
                cell       = [ SX, SY, SZ ]
                ) ]

    # Compute moments
    M = [ ]
    for i in range(8):
        M += [ get_moments(a[i].get_positions(), a[i].get_charges(), L_MAX, r0) ]

    # Construct a composite atoms object
    # and compute the corresponding multipole
    # expansion
    b = ase.Atoms(cell=[ 2*SX, 2*SY, 2*SZ ])
    Mc0_l, Mc_L = zero_moments(L_MAX)
    for i, ( ca, ( M0_l, M_L ) ) in enumerate(zip(a, M)):
        x = i % 2
        y = (i/2) % 2
        z = (i/4) % 2

        dr = np.array([x*SX, y*SY, z*SZ])

        ca.translate(dr)
        b += ca

        dr = np.array([(2*x-1)*SX/2, (2*y-1)*SY/2, (2*z-1)*SZ/2])

        multipole_to_multipole(dr, L_MAX, M0_l, M_L, Mc0_l, Mc_L)

    # Compute the multipole moment directly
    Md0_l, Md_L = get_moments(b.get_positions(), b.get_charges(), L_MAX, r0c)

    assert np.max(np.abs(Mc0_l-Md0_l)) < THRES_MOM and np.max(np.abs(Mc_L-Md_L)) < THRES_MOM

    # Now that we have verified that the moments are okay,
    # lets compute the field somewhere randomly.
    x, y, z = np.random.random_integers(-3,3,3)
    while abs(x) != 3 and abs(y) != 3 and abs(z) != 3:
        x, y, z = np.random.random_integers(-3,3,3)

    r0tar = ( np.array([x, y, z]) + np.random.random(3) ) * np.array([2*SX, 2*SY, 2*SZ])

    L0_l, L_L = multipole_to_local(r0tar - r0c, L_MAX, Mc0_l, Mc_L)

    phi = 0.0
    E   = np.zeros(3)
    for i in b:
        dr = i.get_position() - r0tar
        phi += i.get_charge()/sqrt(np.dot(dr, dr))
        E   -= i.get_charge()*dr/(np.dot(dr, dr)**(3./2))

    assert abs(phi - L0_l[0]) < THRES_PHI
    assert np.max(np.abs(E - np.array([ -L_L[0].real, -L_L[0].imag, L0_l[1] ]))) < THRES_PHI

    # Shift the expansion somewhere else
    r0tar2 = ( np.array([x, y, z]) + np.random.random(3) ) * np.array([2*SX, 2*SY, 2*SZ])

    L0_l2, L_L2 = local_to_local(r0tar2 - r0tar, L_MAX, L0_l, L_L, L_MAX)

    phi = 0.0
    E   = np.zeros(3)
    for i in b:
        dr = i.get_position() - r0tar2
        phi += i.get_charge()/sqrt(np.dot(dr, dr))
        E   -= i.get_charge()*dr/(np.dot(dr, dr)**(3./2))

    assert abs(phi - L0_l2[0]) < THRES_PHI2
    assert np.max(np.abs(E - np.array([ -L_L2[0].real, -L_L2[0].imag, L0_l2[1] ]))) < THRES_PHI2

    # Compute the multipole moment directly
    Md0_l, Md_L = get_moments(b.get_positions(), b.get_charges(), L_MAX, r0c)

    # Now rotate the atoms, recompute the multipole expansion and compare the result
    # to the rotated expansion
    a1 = 2*pi*random.random()
    a2 = 2*pi*random.random()
    a3 = 2*pi*random.random()
    Rz1 = np.array( [ [ cos(a1), -sin(a1), 0 ],
                      [ sin(a1),  cos(a1), 0 ],
                      [       0,        0, 1 ] ] )
    Ry  = np.array( [ [ cos(a2), 0, -sin(a2) ],
                      [       0, 1,        0 ],
                      [ sin(a2), 0,  cos(a2) ] ] )
    Rz2 = np.array( [ [ cos(a3), -sin(a3), 0 ],
                      [ sin(a3),  cos(a3), 0 ],
                      [       0,        0, 1 ] ] )

    R = np.dot(np.dot(Rz1, Ry), Rz2)

    for i in b:
        i.set_position(r0c + np.dot(R, i.get_position() - r0c))

    Me0_l, Me_L = get_moments(b.get_positions(), b.get_charges(), L_MAX, r0c)
    Mf0_l, Mf_L = transform_multipole(R, L_MAX, Md0_l, Md_L)

    assert np.max(np.abs(Me0_l-Mf0_l)) < THRES_ROT and np.max(np.abs(Me_L-Mf_L)) < THRES_ROT


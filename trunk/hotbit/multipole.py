"""
Multipole expansion module

Naming conventions for arrays:

index  description
-----  -----------
l      angular momentum number
L      l and m for l != 0
"""

# Copyright (C) 2010 NSC Jyvaskyla, Fh-IWM
# Please see the accompanying LICENSE file for further information.

import numpy as np

from _hotbit import solid_harmonic_R, multipole_to_multipole
from _hotbit import multipole_to_local, local_to_local, transform_multipole


def lm2index(l, m):
    """
    Return unique index for the off-diagonal elements
    of the multipole moments.

    The mapping is l,m -> l*(l-1)/2 + m - 1
 
    Parameters:
    -----------
    l, m
    """

    return l*(l-1)/2 + m - 1


def zero_moments(l_max):
    """
    Return array with zero moments.
    """

    return np.zeros([l_max+1], dtype=float), \
        np.zeros([lm2index(l_max,l_max)+1], dtype=complex)


def get_moments(r, q, l_max, r0):
    """
    Compute the multipole moments of a set of atoms

    Parameters:
    r:      Positions
    q:      Charges
    l_max:  Maximum angular momentum number for the expansion
    r0:     Expansion origin
    """

    # This could be easily moved to C if it turns out to be slow
    M0_l, M_L  = zero_moments(l_max)
    for r, q in zip(r, q):
        cR0_l, cR_L = solid_harmonic_R(r-r0, l_max)
        M0_l += q*cR0_l
        M_L  += q*cR_L

    return M0_l, M_L.conj()


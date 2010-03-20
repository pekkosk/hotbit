"""
Multipole expansion module

Naming conventions for arrays:

index  description
-----  -----------
l      angular momentum number
L      l and m for l != 0
"""

import numpy as np

from _hotbit import solid_harmonic_R


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

    return np.zeros([l_max+1], dtype=float), np.zeros([lm2index(l_max,l_max)+1], dtype=complex)


def get_moments(a, l_max, r0=None):
    """
    Compute the multipole moments of a set of atoms

    Parameters:
    a:      ASE Atoms object
    l_max:  Maximum angular momentum number for the expansion
    r0:     Expansion origin, if omitted the center of mass will be used
    """

    # This could be easily moved to C if it turns out to be slow

    if r0 is None:
        # FIXME: Should we rather use the position
        # where the number of non-zero multipole
        # moments is minimal?
        r0 = a.get_center_of_mass()

    M0_l = np.zeros([l_max+1], dtype=float)
    M_L  = np.zeros([lm2index(l_max,l_max)+1], dtype=complex)
    for r, q in zip(a.get_positions(), a.get_charges()):
        cR0_l, cR_L = solid_harmonic_R(r-r0, l_max)
        M0_l += q*cR0_l
        M_L  += q*cR_L

    return M0_l, M_L.conj()


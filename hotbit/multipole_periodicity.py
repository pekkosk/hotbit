"""
This module contains the MultipolePeriodicity class, which computes the
electrostatic potential and field for a set of point charges under the
respective symmetry operations.

The potential and field returned here does not include the contribution
of the shape (Gaussian/Slater) of the charges, which is short ranged and
can be easily added later.
"""

import numpy as np

from multipole import get_moments, zero_moments
from multipole import multipole_to_multipole, multipole_to_local
from multipole import local_to_local, transform_multipole

from box.timing import Timer


class MultipolePeriodicity:
    def __init__(self, l_max=8, n=3, k=5, timer=None):
        """
        Instantiate a new MultipolePeriodicity object which computes
        the electrostatic interaction by direct summation using a
        telescoped multipole expansion.

        Parameters:
        -----------
        l_max:   Order of the expansion (maximum angular momentum,
                 typically 5 to 8)
        n:       Number of cells to combine during each telescoping step
        k:       Summation cutoff. The interaction range will be n**k 
                 (number of cells, typically k = 5).
        """
        if l_max < 1:
            raise ValueError("l_max must be >= 1.")
        if n % 2 != 1 or n < 3:
            raise ValueError("k must be >= 3 and odd.")
        if k < 1:
            raise ValueError("k must be >= 1.")

        self.l_max  = l_max
        self.n      = (n-1)/2
        self.k      = k

        if timer is None:
            self.timer  = Timer('MultipolePeriodicity')
        else:
            self.timer  = timer


    def update(self, a, q):
        """
        Compute multipoles.

        Parameters:
        -----------
        r:   Hotbit Atoms object, or atoms object that implements the transform
             and rotation interface.
        q:   Charges
        """
        self.timer.start('multipole_to_multipole')

        r  = a.get_positions()

        # FIXME!!! Center of gravity okay?
        # Probably not, needs to center for
        # rotation operations.
        self.r0  = np.sum(r, 0)/r.shape[0]

        T0_l, T_L = get_moments(r, q, self.l_max, self.r0)

        self.M = [ ( T0_l.copy(), T_L.copy() ) ]

        s1, s2, s3  = a.get_symmetry_operation_ranges()

        # Compute telescoped multipoles
        level = 1
        for k in range(self.k-1):
            M0_l  = T0_l.copy()
            M_L   = T_L.copy()

            for x1 in range(-self.n, self.n+1):
                for x2 in range(-self.n, self.n+1):
                    for x3 in range(-self.n, self.n+1):
                        # Loop over all symmetry operations and compute
                        # telescoped multipoles
                        # FIXME!!! Currently only support continuous symmetries,
                        # think about discrete/recurrent ones.
                        # FIXME!!! Currently assumes periodicity in all three
                        # spatial directions.

                        # The origin is already okay, skip it
                        if x1 != 0 or x2 != 0 or x3 != 0:
                            r = a.transform([0.0,0.0,0.0], [x1,x2,x3])*level
                            multipole_to_multipole(r, self.l_max,
                                                   M0_l, M_L, T0_l, T_L)

            self.M += [ ( T0_l.copy(), T_L.copy() ) ]

            level *= 2*self.n+1

        self.timer.stop('multipole_to_multipole')

        ###

        self.timer.start('multipole_to_local')

        # Compute the local expansion from telescoped multipoles
        L0_l, L_L  = zero_moments(self.l_max)
        n_max      = self.n + 2*self.n+1
        Mi         = len(self.M)
        for k in range(self.k-1):
            level /= 2*self.n+1
            Mi    -= 1

            M0_l, M_L  = self.M[Mi]

            for x1 in range(-n_max, n_max+1):
                for x2 in range(-n_max, n_max+1):
                    for x3 in range(-n_max, n_max+1):
                        # Loop over all symmetry operations and compute the
                        # local expansion from the telescoped multipoles

                        # No local expansion in the inner region
                        if abs(x1) > self.n and abs(x2) > self.n and abs(x3) > self.n:
                            r = a.transform([0.0,0.0,0.0], [x1,x2,x3])*level
                            multipole_to_local(-r, self.l_max, M0_l, M_L, L0_l, L_L)

        self.L = ( L0_l, L_L )

        self.timer.stop('multipole_to_local')

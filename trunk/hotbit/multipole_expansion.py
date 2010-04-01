"""
This module contains the MultipoleExpansion class, which computes the
electrostatic potential and field for a set of point charges under the
respective symmetry operations.

The potential and field returned here does not include the contribution
of the shape (Gaussian/Slater) of the charges, which is short ranged and
can be easily added later.
"""

# Copyright (C) 2010 NSC Jyvaskyla, Fh-IWM
# Please see the accompanying LICENSE file for further information.

from math import pi, sqrt

import numpy as np

from multipole import get_moments, zero_moments
from multipole import multipole_to_multipole, multipole_to_local
from multipole import local_to_local, transform_multipole

from box.timing import Timer


def n_from_ranges(s, n):
    r = [ ]

    for i1, i2 in s:
        if i1 == -np.Inf:
            j1 = -n
        else:
            j1 = int(i1)
        if i2 == np.Inf:
            j2 = n+1
        else:
            j2 = int(i2)+1
        r += [ ( j1, j2 ) ]

    return r


class MultipoleExpansion:
    def __init__(self, l_max=8, n=3, k=5, timer=None):
        """
        Instantiate a new MultipoleExpansion object which computes
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
            raise ValueError("n must be >= 3 and odd.")
        if k < 1:
            raise ValueError("k must be >= 1.")

        self.l_max  = l_max
        self.n      = (n-1)/2
        self.k      = k

        if timer is None:
            self.timer  = Timer('MultipoleExpansion')
        else:
            self.timer  = timer


    def update(self, a, q):
        """
        Compute multipoles, do the transformations, and compute the
        electrostatic potential and field on each atom in a.

        Parameters:
        -----------
        a:   Hotbit Atoms object, or atoms object that implements the transform
             and rotation interface.
        q:   Charges
        """
        self.timer.start('multipole_to_multipole')

        nat  = len(a)
        r    = a.get_positions()

        # FIXME!!! Center of box okay?
        # Probably not, needs to center for
        # rotation operations.
        cell = a.get_cell()
        self.r0  = ( cell[0, :] + cell[1, :] + cell[2, :] )/2

        T0_l, T_L = get_moments(r, q, self.l_max, self.r0)

        self.M = [ ( T0_l.copy(), T_L.copy() ) ]

        sym_ranges  = a.get_symmetry_operation_ranges()
        n1, n2, n3  = n_from_ranges(sym_ranges, self.n)

        # Compute telescoped multipoles
        level = 1
        for k in range(self.k-2):
            M0_l  = T0_l.copy()
            M_L   = T_L.copy()

            for x1 in range(*n1):
                for x2 in range(*n2):
                    for x3 in range(*n3):
                        # Loop over all symmetry operations and compute
                        # telescoped multipoles
                        # FIXME!!! Currently only supports continuous symmetries,
                        # think about discrete/recurrent ones.
                        # FIXME!!! No rotations yet

                        # The origin is already okay, skip it
                        if x1 != 0 or x2 != 0 or x3 != 0:
                            r1 = a.transform(self.r0,
                                             [x1*level, x2*level, x3*level])
                            multipole_to_multipole(r1-self.r0, self.l_max,
                                                   M0_l, M_L, T0_l, T_L)

            self.M += [ ( T0_l.copy(), T_L.copy() ) ]

            level *= 2*self.n+1

        self.timer.stop('multipole_to_multipole')

        ###

        self.timer.start('multipole_to_local')

        # Compute the local expansion from telescoped multipoles
        L0_l, L_L      = zero_moments(self.l_max)
        nm1, nm2, nm3  = n_from_ranges(sym_ranges, self.n + 2*self.n + 1)
        Mi             = len(self.M)-1
        for k in range(self.k-1):
            M0_l, M_L  = self.M[Mi]

            for x1 in range(*nm1):
                for x2 in range(*nm2):
                    for x3 in range(*nm3):
                        # Loop over all symmetry operations and compute the
                        # local expansion from the telescoped multipoles

                        # No local expansion in the inner region
                        if abs(x1) > self.n or abs(x2) > self.n or \
                                abs(x3) > self.n:
                            r1 = a.transform(self.r0,
                                             [x1*level, x2*level, x3*level])
                            multipole_to_local(r1-self.r0, self.l_max,
                                               M0_l, M_L, L0_l, L_L)

            level /= 2*self.n+1
            Mi    -= 1

        self.L = ( L0_l, L_L )

        self.timer.stop('multipole_to_local')

        ###

        self.phi_a  = np.zeros(nat, dtype=float)
        self.E_av    = np.zeros([nat, 3], dtype=float)
    
        ###

        self.timer.start('local_to_local')

        for i in a:
            loc0_l, loc_L          = local_to_local(i.get_position()-self.r0,
                                                   self.l_max, L0_l, L_L, 1)
            self.phi_a[i.index]    = loc0_l[0]
            self.E_av[i.index, :]  = [ -loc_L[0].real,
                                       -loc_L[0].imag,
                                        loc0_l[1] ]

        self.timer.stop('local_to_local')

        ###

        self.timer.start('near_field')

        # Contributing of neighboring boxes
        for x1 in range(*n1):
            for x2 in range(*n2):
                for x3 in range(*n3):
                    # self-interaction needs to be treated separately
                    if x1 != 0 or x2 != 0 or x3 != 0:
                        # construct a matrix with distances
                        r1      = a.transform(self.r0, [x1,x2,x3])

                        dr      = r.reshape(nat, 1, 3) - \
                            (r1+r-self.r0).reshape(1, nat, 3)
                        abs_dr  = np.sqrt(np.sum(dr*dr, axis=2))
                        phi     = q/abs_dr
                        E       = q.reshape(1, nat, 1)*dr/ \
                            (abs_dr**3).reshape(nat, nat, 1)

                        self.phi_a += np.sum(phi, axis=1)
                        self.E_av  += np.sum(E, axis=1)

        # Self-contribution
        dr        = r.reshape(nat, 1, 3) - r.reshape(1, nat, 3)
        abs_dr    = np.sqrt(np.sum(dr*dr, axis=2))

        # Avoid divide by zero
        abs_dr[np.diag_indices_from(abs_dr)]  = 1.0

        phi       = q/abs_dr
        E         = q.reshape(1, nat, 1)*dr/(abs_dr**3).reshape(nat, nat, 1)

        phi[np.diag_indices_from(phi)]   = 0.0
        E[np.diag_indices_from(phi), :]  = 0.0

        self.phi_a += np.sum(phi, axis=1)
        self.E_av  += np.sum(E, axis=1)

        # Dipole correction for 3D sum
        s1, s2, s3 = sym_ranges
        if s1[1] == np.Inf and s2[1] == np.Inf and s3[1] == np.Inf:
            Ml0, Mlm = self.M[0]

            dip  = np.array([-2*Mlm[0].real, 2*Mlm[0].imag, Ml0[1]])
            dip *= 4*pi/(3*a.get_volume())

            self.phi_a -= np.dot(r-self.r0, dip)
            self.E_av  += dip

        self.timer.stop('near_field')


    def get_moments(self):
        """
        Return the multipole moments.
        """
        return self.M

    
    def get_local_expansion(self):
        """
        Return the local expansion of the potential.
        """
        return self.L


    def get_potential(self):
        """
        Return the electrostatic potential for each atom.
        """
        return self.phi_a


    def get_field(self):
        """
        Return the electrostatic field for each atom.
        """
        return self.E_av


    def get_potential_and_field(self):
        """
        Return the both, the electrostatic potential and the field for each atom.
        """
        return self.phi_a, self.E_av

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
from scipy.special import erfc

from hotbit.neighbors import get_neighbors
from box.timing import Timer

from hotbit.coulomb.baseclass import Coulomb

# diag_indices_from was introduced in numpy 1.4.0
if hasattr(np, 'diag_indices_from'):
    diag_indices_from = np.diag_indices_from
else:
    def diag_indices_from(m):
        i = [ ] 
        for n in m.shape:
            i += [ np.arange(n, dtype=int) ]
        return tuple(i)


class DirectCoulomb(Coulomb):
    def __init__(self, cutoff=None, timer=None):
        """
        Instantiate a new DirectCoulomb object which computes the electrostatic
        interaction by direct summation.

        Parameters:
        -----------
        cutoff:   If not None, the Coulomb interaction will be smoothly forced
                  to zero at this distance by multiplication with erfc(r/cutoff)
        """
        if timer is None:
            self.timer  = Timer('DirectCoulomb')
        else:
            self.timer  = timer

        self.cutoff = cutoff

        # Last positions
        self.r_av  = None
        # Last charges
        self.q_a   = None


    def update(self, a, q=None):
        if q is None:
            q = a.get_charges()

        r = a.get_positions()
        # FIXME!!! Check for change in cell, symmetries
        if self.r_av is None or self.q_a is None:
           self._update(a, q)
        elif np.any(r != self.r_av) or np.any(q != self.q_a):
            self._update(a, q)


    def _update(self, a, q):
        """
        Compute the electrostatic potential and field on each atom in a.

        Parameters:
        -----------
        a:   Hotbit Atoms object, or atoms object that implements the transform
             and rotation interface.
        q:   Charges
        """
        self.timer.start('direct_coulomb')

        self.a     = a
        self.r_av  = a.get_positions().copy()
        self.q_a   = q.copy()

        nat  = len(a)

        il, jl, dl, nl = get_neighbors(a, self.cutoff)

        if il is not None:
            if self.cutoff is None:
                phi  = q[jl]/dl
                dl **= 2
                E    = q[jl].reshape(-1, 1)*nl/dl.reshape(-1, 1)
            else:
                f    = erfc(dl/self.cutoff)
                df   = 2/sqrt(pi)*np.exp(-(dl/self.cutoff)**2)/self.cutoff
                phi  = q[jl]*f/dl
                E    = q[jl]*(df + f/dl)/dl
                E    = E.reshape(-1, 1)*nl

        self.phi_a  = np.zeros(nat, dtype=float)
        self.E_av   = np.zeros([nat, 3], dtype=float)

        if il is not None:
            # FIXME!!! Is there some fast numpy magic to compute this?
            for i in xrange(nat):
                self.phi_a[i]    = phi[il == i].sum()
                self.E_av[i, :]  = E[il == i].sum(axis=0)

        self.timer.stop('direct_coulomb')


    def get_potential(self, a=None):
        """
        Return the electrostatic potential for each atom.
        """
        if a is not None:
            self.update(a)

        return self.phi_a


    def get_field(self, a=None):
        """
        Return the electrostatic field for each atom.
        """
        if a is not None:
            self.update(a)

        return self.E_av


    def get_potential_and_field(self, a=None):
        """
        Return the both, the electrostatic potential and the field for each
        atom.
        """
        if a is not None:
            self.update(a)

        return self.phi_a, self.E_av


    def get_gamma(self, a=None):
        """
        Return the gamma correlation matrix, i.e. phi(i) = gamma(i, j)*q(j)
        """
        if a is not None:
            self.update(a)

        self.timer.start('get_gamma')

        nat = len(self.a)

        il, jl, dl, nl = get_neighbors(self.a, self.cutoff)

        if il is None:
            G = None
        else:
            G = np.zeros([nat, nat], dtype=float)
            if self.cutoff is None:
                for i, j, d in zip(il, jl, dl):
                    G[i, j] += 1.0/d
            else:
                for i, j, d in zip(il, jl, dl):
                    G[i, j] += 1.0*erfc(d/self.cutoff)/d

        self.timer.stop('get_gamma')

        return G


### For use as a standalone calculator

    def get_potential_energy(self, a=None):
        """
        Return the Coulomb energy.
        """
        if a is not None:
            self.update(a)

        return np.sum(self.q_a*self.phi_a)/2

    def get_forces(self, a=None):
        """
        Return forces
        """
        if a is not None:
            self.update(a)

        return self.q_a.reshape(-1, 1)*self.E_av

"""
Ewald summation.

The potential and field returned here does not include the contribution
of the shape (Gaussian/Slater) of the charges, which is short ranged and
can be easily added later.
"""

# Copyright (C) 2010 NSC Jyvaskyla, Fh-IWM
# Please see the accompanying LICENSE file for further information.

from math import log, pi, sqrt

import numpy as np
# FIXME!!! Requires scipy, contribute erfc to numpy
from scipy.special import erfc

from ase.units import Hartree, Bohr

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


# FIXME!!! No field yet in Ewald sum
class EwaldSum(Coulomb):
    def __init__(self, accuracy_goal, weight, timer=None):
        self.accuracy_goal  = accuracy_goal
        self.weight         = weight

        if timer is None:
            self.timer  = Timer('EwaldSum')
        else:
            self.timer  = timer


    def update(self, a, q):
        """
        Compute the electrostatic potential.

        Parameters:
        -----------
        a:   Hotbit Atoms object, or atoms object that implements the transform
             and rotation interface.
        q:   Charges
        """
        self.alpha       = (self.weight*pi**3*len(a)/a.get_volume())**(1./3)
        self.sqrt_alpha  = sqrt(self.alpha)

        self.G_cutoff    = 2*sqrt(log(10.0)*self.accuracy_goal*self.alpha)
        self.r_cutoff    = sqrt(log(10.0)*self.accuracy_goal/self.alpha)

        cell_cv      = a.get_cell()
        rec_cell_vc  = np.linalg.inv(cell_cv)
        r_av         = a.get_positions()

        self.timer.start('reciprocal sum')

        # Reciprocal sum
        lx, ly, lz   = np.sqrt(np.sum(rec_cell_vc**2, axis=0))

        maxGx  = int(self.G_cutoff/(2*pi*lx))+1
        maxGy  = int(self.G_cutoff/(2*pi*ly))+1
        maxGz  = int(self.G_cutoff/(2*pi*lz))+1

        Gx  = 2*pi * np.arange(-maxGx, maxGx+1).reshape(-1,  1,  1,  1)
        Gy  = 2*pi * np.arange(-maxGy, maxGy+1).reshape( 1, -1,  1,  1)
        Gz  = 2*pi * np.arange(-maxGz, maxGz+1).reshape( 1,  1, -1,  1)

        G   = Gx*np.array([1,0,0])+Gy*np.array([0,1,0])+Gz*np.array([0,0,1])
        G   = np.dot(G, rec_cell_vc)

        si  = np.sum( np.sin(np.tensordot(G, r_av, axes=(3,1)))*q, axis=3)
        co  = np.sum( np.cos(np.tensordot(G, r_av, axes=(3,1)))*q, axis=3)

        G_sq   = np.sum( G*G, axis=3 )

        rec_G_sq = np.zeros_like(G_sq)
        rec_G_sq[G_sq>0.0] = 1.0/G_sq[G_sq>0.0]
        rec_G_sq[maxGx, maxGy, maxGz] = 0.0

        phase  = np.tensordot(G, r_av, axes=(3, 1))

        si.shape = ( 2*maxGx+1, 2*maxGy+1, 2*maxGz+1, 1 )
        co.shape = ( 2*maxGx+1, 2*maxGy+1, 2*maxGz+1, 1 )

        self.phi_a  = np.sum( np.sum( np.sum(
                    ( np.exp(-G_sq/(4*self.alpha))*rec_G_sq 
                      ).reshape(2*maxGx+1, 2*maxGy+1, 2*maxGz+1, 1) *
                    ( si * np.sin(phase) + co * np.cos(phase) ),
                    axis=0 ), axis=0 ), axis=0 )
        self.phi_a *= 4*pi/a.get_volume()

        self.timer.stop('reciprocal sum')

        self.timer.start('real space sum')

        # Real space sum
        lx, ly, lz   = np.sqrt(np.sum(cell_cv**2, axis=1))

        maxrx  = int(self.r_cutoff/lx)+1
        maxry  = int(self.r_cutoff/ly)+1
        maxrz  = int(self.r_cutoff/lz)+1

        nat  = len(a)
        r    = a.get_positions()
        for x in range(-maxrx, maxrx+1):
            for y in range(-maxry, maxry+1):
                for z in range(-maxrz, maxrz+1):
                    if x != 0 or y != 0 or z != 0:
                        r1          = np.dot([x,y,z], cell_cv)

                        dr          = r.reshape(nat, 1, 3) - \
                            (r1+r).reshape(1, nat, 3)
                        abs_dr      = np.sqrt(np.sum(dr*dr, axis=2))

                        phi         = q*erfc(self.sqrt_alpha*abs_dr)/abs_dr

                        self.phi_a += np.sum(phi, axis=1)
                        

        ## Self-contribution
        dr        = r.reshape(nat, 1, 3) - r.reshape(1, nat, 3)
        abs_dr    = np.sqrt(np.sum(dr*dr, axis=2))

        ## Avoid divide by zero
        abs_dr[diag_indices_from(abs_dr)]  = 1.0

        phi       = q*erfc(self.sqrt_alpha*abs_dr)/abs_dr

        phi[diag_indices_from(phi)]   = 0.0

        self.phi_a += np.sum(phi, axis=1)

        self.timer.stop('real space sum')

        # Self energy
        self.phi_a -= 2*q*sqrt(self.alpha/pi)


    def get_potential(self):
        """
        Return the electrostatic potential for each atom.
        """
        return self.phi_a


### For use as a standalone calculator
### Note: These functions assume eV/A units

    def get_potential_energy(self, a, q=None):
        if q is None:
            q = a.get_charges()

        self.update(a, q)
        return Hartree * Bohr * np.sum(q*self.phi_a)/2

"""
Ewald summation. Reference for the multipole summation.
"""

from math import log, pi, sqrt

import numpy as np

from scipy.special import erfc

from ase.units import Hartree, Bohr


class EwaldSum:
    def __init__(self, accuracy_goal, weight):
        self.accuracy_goal  = accuracy_goal
        self.weight         = weight


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

        rec_G_sq = 1.0/G_sq
        rec_G_sq[maxGx, maxGy, maxGz] = 0.0

        phase  = np.tensordot(G, r_av, axes=(3, 1))

        si.shape = ( 2*maxGx+1, 2*maxGy+1, 2*maxGz+1, 1 )
        co.shape = ( 2*maxGx+1, 2*maxGy+1, 2*maxGz+1, 1 )

        self.phi  = np.sum( np.sum( np.sum(
                    ( np.exp(-G_sq/(4*self.alpha))*rec_G_sq 
                      ).reshape(2*maxGx+1, 2*maxGy+1, 2*maxGz+1, 1) *
                    ( si * np.sin(phase) + co * np.cos(phase) ),
                    axis=0 ), axis=0 ), axis=0 )
        self.phi *= 4*pi/a.get_volume()

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
                        r1        = np.dot([x,y,z], cell_cv)

                        dr        = r.reshape(nat, 1, 3) - \
                            (r1+r).reshape(1, nat, 3)
                        abs_dr    = np.sqrt(np.sum(dr*dr, axis=2))

                        phi       = q*erfc(self.sqrt_alpha*abs_dr)/abs_dr

                        self.phi += np.sum(phi, axis=1)
                        

        ## Self-contribution
        dr        = r.reshape(nat, 1, 3) - r.reshape(1, nat, 3)
        abs_dr    = np.sqrt(np.sum(dr*dr, axis=2))

        ## Avoid divide by zero
        abs_dr[np.diag_indices_from(abs_dr)]  = 1.0

        phi       = q*erfc(self.sqrt_alpha*abs_dr)/abs_dr

        phi[np.diag_indices_from(phi)]   = 0.0

        self.phi += np.sum(phi, axis=1)

        # Self energy
        self.phi -= 2*q*sqrt(self.alpha/pi)

        self.phi *= Hartree*Bohr



    def get_potential(self):
        """
        Return the electrostatic potential for each atom.
        """
        return self.phi

"""
Compute electrostatic interaction and construct the 'electrostatic'
matrix h1 containing the shift in on-site energies.

This module only adds the contribution of the charge density. The interaction
of point charges is computed by one of the solvers in the Coulomb submodule.
"""

# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

from math import exp, log, pi, sqrt
from weakref import proxy
from math import sqrt

from ase.units import Hartree, Bohr

import numpy as np
from scipy.special import erf, erfc
norm=np.linalg.norm
dot=np.dot

from neighbors import get_neighbors
from hotbit.coulomb import DirectCoulomb

# diag_indices_from was introduced in numpy 1.4.0
if hasattr(np, 'diag_indices_from'):
    diag_indices_from = np.diag_indices_from
else:
    def diag_indices_from(m):
        i = [ ] 
        for n in m.shape:
            i += [ np.arange(n, dtype=int) ]
        return tuple(i)


# Some constants
log2 = log(2.0)


def get_Gaussian_gamma_correction(a, U, FWHM=None,
                                  cutoff=None, accuracy_goal=12):
    """
    Gaussian charge distribution.
    """
    if FWHM is None:
        FWHM = sqrt(8*log2/pi)/U

    max_FWHM = np.max(FWHM)
    if max_FWHM <= 0.0:
        raise ValueError("Maximum FWHM (%f) smaller than or equal to zero. " %
                         max_FWHM)

    if cutoff is None and not a.is_cluster():
        # Estimate a cutoff from the accuracy goal if the
        # system is periodic and no cutoff was given
        cutoff = sqrt(log(10.0)*accuracy_goal*max_FWHM/(sqrt(4*log2)))

    il, jl, dl, nl = get_neighbors(a, cutoff)

    nat = len(a)
    G   = np.zeros([nat, nat], dtype=float)
    dG  = np.zeros([nat, nat, 3], dtype=float)
    G[diag_indices_from(G)] = U

    if il is not None:
        for i, j, d, n in zip(il, jl, dl/Bohr, nl):
            const        = 2*sqrt( log2/(FWHM[i]**2+FWHM[j]**2) )
            ecr          = -erfc(const*d)
            decr         = 2/sqrt(pi)*exp(-(const*d)**2)*const
            G[i, j]     += ecr/d
            dG[i, j, :] += (ecr/d-decr)*n/d

    return G, dG


def get_Slater_gamma_correction(a, U, FWHM=None,
                                cutoff=None, accuracy_goal=12):
    """
    Slater-type charge distribution.
    See: M. Elstner et al., Phys. Rev. B 58, 7260 (1998)
    """
    min_U = np.min(U)
    if min_U <= 0.0:
        raise ValueError("Minimum U (%f) smaller than or equal to zero. " %
                         min_U)

    if cutoff is None and not a.is_cluster():
        # Estimate a cutoff from the accuracy goal if the
        # system is periodic and no cutoff was given
        cutoff = sqrt(log(10.0)*accuracy_goal/(sqrt(pi/2)*min_U))

    tau = 16*np.asarray(U)/5

    il, jl, dl, nl = get_neighbors(a, cutoff)

    nat = len(a)
    G   = np.zeros([nat, nat], dtype=float)
    dG  = np.zeros([nat, nat, 3], dtype=float)
    G[diag_indices_from(G)] = U

    if il is not None:
        for i, j, d, n in zip(il, jl, dl/Bohr, nl):
            if abs(tau[i] - tau[j]) < 1e-6:
                src = 1.0/(tau[i]+tau[j])
                fac = tau[i]*tau[j]*src
                avg = 1.6*(fac+fac*fac*src)
                fac = avg*d
                fac2 = fac*fac
                efac = exp(-fac)/(48*d)
                h = -(48 + 33*fac + fac2*(9+fac))*efac
                G[i, j] += \
                    h
                dG[i, j, :] += \
                    (   h/d \
                      + avg*h \
                      + (33*avg + 18*fac*avg + 3*fac2*avg)*efac )*n
            else:
                fi1  = 1.0/(2*(tau[i]**2-tau[j]**2)**2)
                fj1  = -tau[i]**4*tau[j]*fi1
                fi1 *= -tau[j]**4*tau[i]

                fi2  = 1.0/((tau[i]**2-tau[j]**2)**3)
                fj2  = -(tau[i]**6-3*tau[i]**4*tau[j]**2)*fi2
                fi2 *=  (tau[j]**6-3*tau[j]**4*tau[i]**2)

                expi = exp(-tau[i]*d)
                expj = exp(-tau[j]*d)

                G[i, j] += \
                      expi*(fi1+fi2/d) \
                    + expj*(fj1+fj2/d)
                dG[i, j, :] += \
                    (   expi*(tau[i]*(fi1+fi2/d) + fi2/(d**2)) \
                      + expj*(tau[j]*(fj1+fj2/d) + fj2/(d**2)) )*n

    return G, dG


_gamma_correction_dict = {
    'Gaussian': get_Gaussian_gamma_correction,
    'Slater': get_Slater_gamma_correction
    }


class Electrostatics:
    def __init__(self, calc, charge_density='Gaussian',
                 solver=None, accuracy_goal=12):
        self.calc=proxy(calc)
        self.norb=calc.el.get_nr_orbitals()
        self.SCC=calc.get('SCC')
        self.N=len(calc.el)
        self.dq=np.zeros((self.N))
        self.G=np.zeros((self.N,self.N))
        self.dG=np.zeros((self.N,self.N,3))

        self.accuracy_goal = accuracy_goal

        if not charge_density in _gamma_correction_dict.keys():
            raise RuntimeError("Unknown charge density type: %s." %
                               charge_density)

        self.gamma_correction = _gamma_correction_dict[charge_density]

        if solver is None:
            self.solver = DirectCoulomb(self.calc.get('gamma_cut'))
        else:
            self.solver = solver


    def set_dq(self,dq):
        """ (Re)setting dq gives new gamma-potential. """
        self.dq=dq
        self.epsilon = np.sum(self.G*self.dq, axis=1)
        if self.solver is not None:
            self.solver.update(self.calc.el.atoms, dq)
            # Unit mess: The Coulomb solver is unit agnostic, but the Elements
            # object returns the distances in Bohr (from get_distances())
            self.epsilon += self.solver.get_potential()*Bohr


    def coulomb_energy(self):
        """ Return Coulomb energy. """
        self.calc.start_timing('ecoul')
        ecoul = 0.0
        if not self.SCC:
            ecoul = 0.0
        else:
            ecoul = 0.5 * dot(self.dq,self.epsilon)
        self.calc.stop_timing('ecoul')
        return ecoul


    def gamma_forces(self):
        """ Return forces due to electrostatic interactions. """
        if not self.SCC:
            return np.zeros((self.N,3))
        else:
            self.calc.start_timing('f_es')
            if self.solver is None:
                E = np.zeros((self.N,3))
            else:
                # Unit mess: The Coulomb solver is unit agnostic, but the
                # Elements object returns the distances in Bohr
                # (from get_distances())
                E = self.solver.get_field()*Bohr**2
            depsilon  = np.sum(self.dq.reshape(1, -1, 1)*self.dG, axis=1)
            f         = self.dq.reshape(-1, 1) * ( depsilon + E )
            self.calc.stop_timing('f_es')
            return f


    def construct_h1(self,dq=None):
        """ Make the electrostatic part of the Hamiltonian. """
        self.calc.start_timing('h1')
        if dq!=None:
            self.set_dq(dq)
                          
        lst = self.calc.el.get_property_lists(['i','o1','no'])
        
        aux = np.zeros((self.norb,self.norb))
        for i,o1i,noi in lst:
            aux[o1i:o1i+noi,:] = self.epsilon[i]
        h1 = 0.5 * ( aux+aux.transpose() )
                
        # external electrostatics
        if np.any(abs(self.ext)>1E-12):
            aux = np.zeros((self.norb,self.norb))
            for i,o1i,noi in lst:
                aux[o1i:o1i+noi,:] = self.ext[i]
            h1 = h1 + 0.5 * (-1) * (aux+aux.transpose())
                
        self.calc.stop_timing('h1')
        self.h1=h1
        return self.h1


    def get_h1(self):
        """ Get the current electrostatic Hamiltonian. """
        return self.h1


    def construct_Gamma_matrix(self, a):
        '''
        Construct the G-matrix and its derivative.
        
        Done once for each geometry; 
        G_ij = sum_n gamma_ij(Rijn)
        dG_ij = sum_n gamma'_ij(Rijn) hat(Rijn) 
        '''
        self.calc.start_timing('gamma matrix')

        U = [ self.calc.el.get_element(i).get_U()
              for i in range(len(a)) ]
        FWHM = [ self.calc.el.get_element(i).get_FWHM()
                 for i in range(len(a)) ]

        G, dG = self.gamma_correction(a, U,
                                      FWHM = FWHM,
                                      cutoff = self.calc.get('gamma_cut'),
                                      accuracy_goal = self.accuracy_goal)

        self.G, self.dG  = G, dG

        self.ext = np.array( [self.calc.env.phi(i) for i in range(self.N)] )
        self.calc.stop_timing('gamma matrix')


    def get_gamma(self):
        return self.G + self.solver.get_gamma()*Bohr


### For use as a standalone calculator, return eV/A units

    def get_potential_energy(self, a):
        self.calc.el.update_geometry(a)
        self.construct_Gamma_matrix(a)
        self.set_dq(a.get_charges())
        return self.coulomb_energy()*Hartree

    def get_forces(self, a):
        self.calc.el.update_geometry(a)
        self.construct_Gamma_matrix(a)
        self.set_dq(a.get_charges())
        return self.gamma_forces()*Hartree/Bohr
        

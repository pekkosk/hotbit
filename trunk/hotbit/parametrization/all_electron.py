# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""
Atomic Density Functional Theory
"""

from math import pi, sqrt, log
import pickle
import sys

import numpy as npy
from ase.data import atomic_names

from gpaw.atom.configurations import configurations
from gpaw.xc_functional import XCRadialGrid, XCFunctional
from gpaw.utilities import hartree, devnull
from gpaw.exx import atomic_exact_exchange

# fine-structure constant
alpha = 1 / 137.036

class RadialGridDescriptor:
    """Descriptor-class for radial grid."""
    def __init__(self, r_g, dr_g):
        """Construct `RadialGridDescriptor`.

        The one-dimensional array ``r_g`` gives the radii of the grid
        points according to some possibly non-linear function:
        ``r_g[g]`` = *f(g)*.  The array ``dr_g[g]`` = *f'(g)* is used
        for forming derivatives."""
        
        self.r_g = r_g
        self.dr_g = dr_g
        self.dv_g = 4 * pi * r_g**2 * dr_g

    def derivative(self, n_g, dndr_g):
        """Finite-difference derivative of radial function."""
        dndr_g[0] = n_g[1] - n_g[0]
        dndr_g[1:-1] = 0.5 * (n_g[2:] - n_g[:-2])
        dndr_g[-1] = n_g[-1] - n_g[-2]
        dndr_g /= self.dr_g

    def derivative2(self, a_g, b_g):
        """Finite-difference derivative of radial function.

        For an infinitely dense grid, this method would be identical
        to the `derivative` method."""
        
        c_g = a_g / self.dr_g
        b_g[0] = 0.5 * c_g[1] + c_g[0]
        b_g[1] = 0.5 * c_g[2] - c_g[0]
        b_g[1:-1] = 0.5 * (c_g[2:] - c_g[:-2])
        b_g[-2] = c_g[-1] - 0.5 * c_g[-3]
        b_g[-1] = -c_g[-1] - 0.5 * c_g[-2]


class AllElectron:
    """Object for doing an atomic DFT calculation."""

    def __init__(self, symbol, xcname='LDA', scalarrel=False,
                 corehole=None, configuration=None, nofiles=True,
                 txt='-'):
        """Do an atomic DFT calculation.

        Example::

          a = AllElectron('Fe')
          a.run()
        """
        
        if txt is None:
            txt = devnull
        elif txt == '-':
            txt = sys.stdout
        elif isinstance(txt, str):
            txt = open(txt, 'w')
        self.txt = txt

        self.symbol = symbol
        self.xcname = xcname
        self.scalarrel = scalarrel
        self.nofiles = nofiles

        # Get reference state:
        self.Z, nlfe_j = configurations[symbol]

        # Collect principal quantum numbers, angular momentum quantum
        # numbers, occupation numbers and eigenvalues (j is a combined
        # index for n and l):
        self.n_j = [n for n, l, f, e in nlfe_j]
        self.l_j = [l for n, l, f, e in nlfe_j]
        self.f_j = [f for n, l, f, e in nlfe_j]
        self.e_j = [e for n, l, f, e in nlfe_j]

        if configuration is not None:
            j = 0
            for conf in configuration.split(','):
                if conf[0].isdigit():
                    n = int(conf[0])
                    l = 'spdf'.find(conf[1])
                    if len(conf) == 2:
                        f = 1.0
                    else:
                        f = float(conf[2:])
                    assert n == self.n_j[j]
                    assert l == self.l_j[j]
                    self.f_j[j] = f
                    j += 1
                else:
                    j += {'He': 1,
                          'Ne': 3,
                          'Ar': 5,
                          'Kr': 8,
                          'Xe': 11}[conf]

        maxnodes = max([n - l - 1 for n, l in zip(self.n_j, self.l_j)])
        self.N = (maxnodes + 1) * 150
        self.beta = 0.4

        t = self.text
        t()
        if scalarrel:
            t('Scalar-relativistic atomic ', end='')
        else:
            t('Atomic ', end='')
        t('%s calculation for %s (%s, Z=%d)' % (xcname, symbol,
                                                atomic_names[self.Z], self.Z))

        if corehole is not None:
            self.ncorehole, self.lcorehole, self.fcorehole = corehole

            # Find j for core hole and adjust occupation:
            for j in range(len(self.f_j)):
                if self.n_j[j] == self.ncorehole and self.l_j[j] == self.lcorehole:
                    assert self.f_j[j] == 2 * (2 * self.lcorehole + 1)
                    self.f_j[j] -= self.fcorehole
                    self.jcorehole = j
                    break

            coreholestate='%d%s' % (self.ncorehole, 'spd'[self.lcorehole])
            t('Core hole in %s state (%s occupation: %.1f)' % (
                coreholestate, coreholestate, self.f_j[self.jcorehole]))
        else:
            self.jcorehole = None
            self.fcorehole = 0

    def text(self, *args, **kwargs):
        self.txt.write(kwargs.get('sep', ' ').join([str(arg)
                                                    for arg in args]) +
                       kwargs.get('end', '\n'))

    def intialize_wave_functions(self):
        r = self.r
        dr = self.dr
        # Initialize with Slater function:
        for l, e, u in zip(self.l_j, self.e_j, self.u_j):
            a = sqrt(-2.0 * e)

            # This one: "u[:] = r**(1 + l) * npy.exp(-a * r)" gives
            # OverflowError: math range error XXX
            u[:] = r**(1 + l)
            rmax = 350.0 / a     # numpy!
            gmax = int(rmax * self.N / (self.beta + rmax))
            u[:gmax] *= npy.exp(-a * r[:gmax])
            u[gmax:] = 0.0

            norm = npy.dot(u**2, dr)
            u *= 1.0 / sqrt(norm)

    def run(self):
        #     beta g
        # r = ------, g = 0, 1, ..., N - 1
        #     N - g
        #
        #        rN
        # g = --------
        #     beta + r

        t = self.text
        N = self.N
        beta = self.beta
        t(N, 'radial gridpoints.')
        g = npy.arange(N, dtype=float)
        self.r = beta * g / (N - g)
        self.dr = beta * N / (N - g)**2
        self.rgd = RadialGridDescriptor(self.r, self.dr)
        self.d2gdr2 = -2 * N * beta / (beta + self.r)**3

        # Number of orbitals:
        nj = len(self.n_j)

        # Radial wave functions multiplied by radius:
        self.u_j = npy.zeros((nj, self.N))

        # Effective potential multiplied by radius:
        self.vr = npy.zeros(N)

        # Electron density:
        self.n = npy.zeros(N)

        # Always spinpaired nspins=1
        self.xcfunc = XCFunctional(self.xcname, 1)
        self.xc = XCRadialGrid(self.xcfunc, self.rgd)

        n_j = self.n_j
        l_j = self.l_j
        f_j = self.f_j
        e_j = self.e_j

        Z = self.Z    # nuclear charge
        r = self.r    # radial coordinate
        dr = self.dr  # dr/dg
        n = self.n    # electron density
        vr = self.vr  # effective potential multiplied by r

        vHr = npy.zeros(self.N)
        self.vXC = npy.zeros(self.N)

        try:
            f = open(self.symbol + '.restart', 'r')
        except IOError:
            self.intialize_wave_functions()
            n[:] = self.calculate_density()
        else:
            if not self.xc.is_non_local():
                t('Using old density for initial guess.')
                n[:] = pickle.load(f)
                n *= Z / (npy.dot(n * r**2, dr) * 4 * pi)
            else:
                # Do not start from initial guess when doing
                # non local XC!
                # This is because we need wavefunctions as well
                # on the first iteration.
                self.intialize_wave_functions()
                n[:] = self.calculate_density()

        bar = '|------------------------------------------------|'
        t(bar)
        niter = 0
        qOK = log(1e-10)
        while True:
            # calculate hartree potential
            hartree(0, n * r * dr, self.beta, self.N, vHr)

            # add potential from nuclear point charge (v = -Z / r)
            vHr -= Z

            # calculated exchange correlation potential and energy
            self.vXC[:] = 0.0

            if self.xc.is_non_local():
                Exc = self.xc.get_non_local_energy_and_potential(self.u_j, self.f_j, self.e_j, self.l_j, self.vXC)
            else:
                tau = None
                if self.xc.xcfunc.mgga:
                    tau = self.calculate_kinetic_energy_density()
                Exc = self.xc.get_energy_and_potential(n, self.vXC)

            # calculate new total Kohn-Sham effective potential and
            # admix with old version
            vr[:] = vHr + self.vXC * r
            if niter > 0:
                vr[:] = 0.4 * vr + 0.6 * vrold
            vrold = vr.copy()

            # solve Kohn-Sham equation and determine the density change
            self.solve()
            dn = self.calculate_density() - n
            n += dn

            # estimate error from the square of the density change integrated
            q = log(npy.sum((r * dn)**2))

            # print progress bar
            if niter == 0:
                q0 = q
                b0 = 0
            else:
                b = int((q0 - q) / (q0 - qOK) * 50)
                if b > b0:
                    self.txt.write(bar[b0:min(b, 50)])
                    self.txt.flush()
                    b0 = b

            # check if converged and break loop if so
            if q < qOK:
                self.txt.write(bar[b0:])
                self.txt.flush()
                break

            niter += 1
            if niter > 117:
                raise RuntimeError, 'Did not converge!'

##         print
        tau = self.calculate_kinetic_energy_density()
##         print "Ekin(tau)=",npy.dot(tau *r**2 , dr) * 4*pi
##         self.write(tau,'tau1')
##         tau2 = self.calculate_kinetic_energy_density2()
##         self.write(tau2,'tau2')
##         self.write(tau-tau2,'tau12')
##         print "Ekin(tau2)=",npy.dot(tau2 *r**2 , dr) * 4*pi

        # When iterations are over calculate the correct exchange energy
        #if self.xc.is_non_local():
        #    Exc = atomic_exact_exchange(self)

        t()
        t('Converged in %d iteration%s.' % (niter, 's'[:niter != 1]))

        if not self.nofiles:
            pickle.dump(n, open(self.symbol + '.restart', 'w'))

        Epot = 2 * pi * npy.dot(n * r * (vHr - Z), dr)
        Ekin = -4 * pi * npy.dot(n * vr * r, dr)
        for f, e in zip(f_j, e_j):
            Ekin += f * e


        t()
        t('Energy contributions:')
        t('-------------------------')
        t('Kinetic:   %+13.6f' % Ekin)
        t('XC:        %+13.6f' % Exc)
        t('Potential: %+13.6f' % Epot)
        t('-------------------------')
        t('Total:     %+13.6f' % (Ekin + Exc + Epot))
        self.ETotal = Ekin + Exc + Epot
        t()

        t('state      eigenvalue         ekin         rmax')
        t('-----------------------------------------------')
        for m, l, f, e, u in zip(n_j, l_j, f_j, e_j, self.u_j):
            # Find kinetic energy:
            k = e - npy.sum((npy.where(abs(u) < 1e-160, 0, u)**2 * #XXXNumeric!
                             vr * dr)[1:] / r[1:])

            # Find outermost maximum:
            g = self.N - 4
            while u[g - 1] >= u[g]:
                g -= 1
            x = r[g - 1:g + 2]
            y = u[g - 1:g + 2]
            A = npy.transpose(npy.array([x**i for i in range(3)]))
            c, b, a = npy.linalg.solve(A, y)
            assert a < 0.0
            rmax = -0.5 * b / a

            s = 'spdf'[l]
            t('%d%s^%-4.1f: %12.6f %12.6f %12.3f' % (m, s, f, e, k, rmax))
        t('-----------------------------------------------')
        t('(units: Bohr and Hartree)')

        for m, l, u in zip(n_j, l_j, self.u_j):
            self.write(u, 'ae', n=m, l=l)

        self.write(n, 'n')
        self.write(vr, 'vr')
        self.write(vHr, 'vHr')
        self.write(self.vXC, 'vXC')
        self.write(tau, 'tau')

        self.Ekin = Ekin
        self.Epot = Epot
        self.Exc = Exc

#mathiasl
       # for x in range(npy.size(self.r)):
       #     print self.r[x] , self.u_j[self.jcorehole,x]


    def write(self, array, name=None, n=None, l=None):
        if self.nofiles:
            return

        if name:
            name = self.symbol + '.' + name
        else:
            name = self.symbol

        if l is not None:
            assert n is not None
            if n > 0:
                # Bound state:
                name += '.%d%s' % (n, 'spdf'[l])
            else:
                name += '.x%d%s' % (-n, 'spdf'[l])

        f = open(name, 'w')
        for r, a in zip(self.r, array):
            print >> f, r, a

    def calculate_density(self):
        """Return the electron charge density divided by 4 pi"""
        n = npy.dot(self.f_j,
                    npy.where(abs(self.u_j) < 1e-160, 0,
                              self.u_j)**2) / (4 * pi)
        n[1:] /= self.r[1:]**2
        n[0] = n[1]
        return n

    def calculate_kinetic_energy_density(self):
        """Return the kinetic energy density"""
        return self.radial_kinetic_energy_density(self.f_j,self.l_j,
                                                  self.u_j)

    def radial_kinetic_energy_density(self,f_j,l_j,u_j):
        """Kinetic energy density from a restricted set of wf's
        """
        shape = npy.shape(u_j[0])
        dudr = npy.zeros(shape)
        tau = npy.zeros(shape)
        for f, l, u in zip(f_j,l_j,u_j):
            self.rgd.derivative(u,dudr)
            # contribution from angular derivatives
            if l>0:
                tau += f * l*(l+1) * npy.where(abs(u) < 1e-160, 0, u)**2
            # contribution from radial derivatives
            dudr = u - self.r*dudr
            tau += f * npy.where(abs(dudr) < 1e-160, 0, dudr)**2
        tau[1:] /= self.r[1:]**4
        tau[0] = tau[1]

        return 0.5 * tau / (4 * pi)

    def calculate_kinetic_energy_density2(self):
        """Return the kinetic energy density
        calculation over R(r)=u(r)/r
        slower convergence with # of radial grid points for
        Ekin of H than radial_kinetic_energy_density
        """

        shape = self.u_j.shape[1]
        R = npy.zeros(shape)
        dRdr = npy.zeros(shape)
        tau = npy.zeros(shape)
        for f, l, u in zip(self.f_j,self.l_j,self.u_j):
            R[1:] = u[1:] / self.r[1:]
            if l==0:
                # estimate value at origin by Taylor series to first order
                d1=self.r[1]
                d2=self.r[2]
                R[0] = .5*(R[1]+R[2]+(R[1]-R[2])*(d1+d2)/(d2-d1))
            else:    R[0] = 0
            self.rgd.derivative(R,dRdr)
            # contribution from radial derivatives
            tau += f * npy.where(abs(dRdr) < 1e-160, 0, dRdr)**2
            # contribution from angular derivatives
            if l>0:
                R[1:] = R[1:] / self.r[1:]
                if l==1: R[0] = R[1]
                else:    R[0] = 0
                tau += f * l*(l+1) * npy.where(abs(R) < 1e-160, 0, R)**2

        return 0.5 * tau / (4 * pi)

    def solve(self):
        """Solve the Schrodinger equation

        ::

             2
            d u     1  dv  du   u     l(l + 1)
          - --- - ---- -- (-- - -) + [-------- + 2M(v - e)] u = 0,
              2      2 dr  dr   r         2
            dr    2Mc                    r


        where the relativistic mass::

                   1
          M = 1 - --- (v - e)
                    2
                  2c

        and the fine-structure constant alpha = 1/c = 1/137.036
        is set to zero for non-scalar-relativistic calculations.

        On the logaritmic radial grids defined by::

              beta g
          r = ------,  g = 0, 1, ..., N - 1
              N - g

                 rN
          g = --------, r = [0; oo[
              beta + r

        the Schrodinger equation becomes::

           2
          d u      du
          --- c  + -- c  + u c  = 0
            2  2   dg  1      0
          dg

        with the vectors c0, c1, and c2  defined by::

                 2 dg 2
          c  = -r (--)
           2       dr

                  2         2
                 d g  2    r   dg dv
          c  = - --- r  - ---- -- --
           1       2         2 dr dr
                 dr       2Mc

                                    2    r   dv
          c  = l(l + 1) + 2M(v - e)r  + ---- --
           0                               2 dr
                                        2Mc
        """
        r = self.r
        dr = self.dr
        vr = self.vr

        c2 = -(r / dr)**2
        c10 = -self.d2gdr2 * r**2 # first part of c1 vector

        if self.scalarrel:
            self.r2dvdr = npy.zeros(self.N)
            self.rgd.derivative(vr, self.r2dvdr)
            self.r2dvdr *= r
            self.r2dvdr -= vr
        else:
            self.r2dvdr = None

        # solve for each quantum state separately
        for j, (n, l, e, u) in enumerate(zip(self.n_j, self.l_j,
                                             self.e_j, self.u_j)):
            nodes = n - l - 1 # analytically expected number of nodes
            delta = -0.2 * e
            nn, A = shoot(u, l, vr, e, self.r2dvdr, r, dr, c10, c2,
                          self.scalarrel)
            # adjust eigenenergy until u has the correct number of nodes
            while nn != nodes:
                diff = cmp(nn, nodes)
                while diff == cmp(nn, nodes):
                    e -= diff * delta
                    nn, A = shoot(u, l, vr, e, self.r2dvdr, r, dr, c10, c2,
                                  self.scalarrel)
                delta /= 2

            # adjust eigenenergy until u is smooth at the turning point
            de = 1.0
            while abs(de) > 1e-9:
                norm = npy.dot(npy.where(abs(u) < 1e-160, 0, u)**2, dr)
                u *= 1.0 / sqrt(norm)
                de = 0.5 * A / norm
                x = abs(de / e)
                if x > 0.1:
                    de *= 0.1 / x
                e -= de
                assert e < 0.0
                nn, A = shoot(u, l, vr, e, self.r2dvdr, r, dr, c10, c2,
                              self.scalarrel)
            self.e_j[j] = e
            u *= 1.0 / sqrt(npy.dot(npy.where(abs(u) < 1e-160, 0, u)**2, dr))


    def solve_confined(self, j, rc, vconf=None):
        """Solve the Schroedinger equation in a confinement potential.
        
        Solves the Schroedinger equation like the solve method, but with a
        number of differences.  Before invoking this method, run solve() to
        get initial guesses.

        Parameters:
            j: solves only for the state given by j
            rc: solution cutoff. Solution will be zero outside this.
            vconf: added to the potential (use this as confinement potential)

        Returns: a tuple containing the solution u and its energy e. 

        Unlike the solve method, this method will not alter any attributes of
        this object.
        """
        r = self.r
        dr = self.dr
        vr = self.vr.copy()
        if vconf is not None:
            vr += vconf

        c2 = -(r / dr)**2
        c10 = -self.d2gdr2 * r**2 # first part of c1 vector

        if j is None:
            n,l,e,u=3,2,-0.15,self.u_j[-1].copy()
        else:
            n = self.n_j[j]
            l = self.l_j[j]
            e = self.e_j[j]
            u = self.u_j[j].copy()
            
        nn, A = shoot_confined(u, l, vr, e, self.r2dvdr, r, dr, c10, c2,
                       self.scalarrel, rc=rc, beta=self.beta)
        assert nn == n - l - 1 # run() should have been called already
        
        # adjust eigenenergy until u is smooth at the turning point
        de = 1.0
        while abs(de) > 1e-9:
            norm = npy.dot(npy.where(abs(u) < 1e-160, 0, u)**2, dr)
            u *= 1.0 / sqrt(norm)
            de = 0.5 * A / norm
            x = abs(de / e)
            if x > 0.1:
                de *= 0.1 / x
            e -= de
            assert e < 0.0

            nn, A = shoot_confined(u, l, vr, e, self.r2dvdr, r, dr, c10, c2,
                                   self.scalarrel, rc=rc, beta=self.beta)
        u *= 1.0 / sqrt(npy.dot(npy.where(abs(u) < 1e-160, 0, u)**2, dr))
        return u, e

    def kin(self, l, u, e=None): # XXX move to Generator
        r = self.r[1:]
        dr = self.dr[1:]

        c0 = 0.5 * l * (l + 1) / r**2
        c1 = -0.5 * self.d2gdr2[1:]
        c2 = -0.5 * dr**-2

        if e is not None and self.scalarrel:
            x = 0.5 * alpha**2
            Mr = r * (1.0 + x * e) - x * self.vr[1:]
            c0 += ((Mr - r) * (self.vr[1:] - e * r) +
                   0.5 * x * self.r2dvdr[1:] / Mr) / r**2
            c1 -= 0.5 * x * self.r2dvdr[1:] / (Mr * dr * r)

        fp = c2 + 0.5 * c1
        fm = c2 - 0.5 * c1
        f0 = c0 - 2 * c2
        kr = npy.zeros(self.N)
        kr[1:] = f0 * u[1:] + fm * u[:-1]
        kr[1:-1] += fp[:-1] * u[2:]
        kr[0] = 0.0
        return kr

    def r2g(self, r):
        """Convert radius to index of the radial grid."""
        return int(r * self.N / (self.beta + r))

    def get_confinement_potential(self, alpha, ri, rc):
        """Create a smooth confinement potential.
        
        Returns a (potential) function which is zero inside the radius ri
        and goes to infinity smoothly at rc, after which point it is nan.
        The potential is given by::

                   alpha         /   rc - ri \
          V(r) = --------   exp ( - --------- )   for   ri < r < rc
                  rc - r         \    r - ri /

        """
        i_ri = self.r2g(ri)
        i_rc = self.r2g(rc)
        if self.r[i_rc] == rc:
            # Avoid division by zero in the odd case that rc coincides
            # exactly with a grid point (which actually happens sometimes)
            i_rc -= 1

        potential = npy.zeros(npy.shape(self.r))
        r = self.r[i_ri+1:i_rc+1]
        exponent = - (rc - ri) / (r - ri)
        denom = rc - r
        value = npy.exp(exponent) / denom
        potential[i_ri+1:i_rc+1] = value
        potential[i_rc+1:] = float('nan')

        return alpha * potential

def shoot(u, l, vr, e, r2dvdr, r, dr, c10, c2, scalarrel=False, gmax=None):
    """n, A = shoot(u, l, vr, e, ...)

    For guessed trial eigenenergy e, integrate the radial Schrodinger
    equation::

          2
         d u      du
         --- c  + -- c  + u c  = 0
           2  2   dg  1      0
         dg

               2 dg 2
        c  = -r (--)
         2       dr

                2         2
               d g  2    r   dg dv
        c  = - --- r  - ---- -- --
         1       2         2 dr dr
               dr       2Mc

                                  2    r   dv
        c  = l(l + 1) + 2M(v - e)r  + ---- --
         0                               2 dr
                                      2Mc

    The resulting wavefunction is returned in input vector u.
    The number of nodes of u is returned in attribute n.
    Returned attribute A, is a measure of the size of the derivative
    discontinuity at the classical turning point.
    The trial energy e is correct if A is zero and n is the correct number
    of nodes."""

    if scalarrel:
        x = 0.5 * alpha**2 # x = 1 / (2c^2)
        Mr = r * (1.0 + x * e) - x * vr
    else:
        Mr = r
    c0 = l * (l + 1) + 2 * Mr * (vr - e * r)
    if gmax is None and npy.alltrue(c0 > 0):
        print """
Problem with initial electron density guess!  Try to run the program
with the '-nw' option (non-scalar-relativistic calculation + write
density) and then try again without the '-n' option (this will
generate a good initial guess for the density).
"""
        
        raise SystemExit
    c1 = c10
    if scalarrel:
        c0 += x * r2dvdr / Mr
        c1 = c10 - x * r * r2dvdr / (Mr * dr)

    # vectors needed for numeric integration of diff. equation
    fm = 0.5 * c1 - c2
    fp = 0.5 * c1 + c2
    f0 = c0 - 2 * c2

    if gmax is None:
        # set boundary conditions at r -> oo (u(oo) = 0 is implicit)
        u[-1] = 1.0

        # perform backwards integration from infinity to the turning point
        g = len(u) - 2
        u[-2] = u[-1] * f0[-1] / fm[-1]
        while c0[g] > 0.0: # this defines the classical turning point
            u[g - 1] = (f0[g] * u[g] + fp[g] * u[g + 1]) / fm[g]
            if u[g - 1] < 0.0:
                # There should't be a node here!  Use a more negative
                # eigenvalue:
                print '!!!!!!',
                return 100, None
            if u[g - 1] > 1e100:
                u *= 1e-100
            g -= 1

        # stored values of the wavefunction and the first derivative
        # at the turning point
        gtp = g + 1
        utp = u[gtp]
        dudrplus = 0.5 * (u[gtp + 1] - u[gtp - 1]) / dr[gtp]
    else:
        gtp = gmax

    # set boundary conditions at r -> 0
    u[0] = 0.0
    u[1] = 1.0

    # perform forward integration from zero to the turning point
    g = 1
    nodes = 0
    while g <= gtp: # integrate one step further than gtp
                    # (such that dudr is defined in gtp)
        u[g + 1] = (fm[g] * u[g - 1] - f0[g] * u[g]) / fp[g]
        if u[g + 1] * u[g] < 0:
            nodes += 1
        g += 1
    if gmax is not None:
        return

    # scale first part of wavefunction, such that it is continuous at gtp
    u[:gtp + 2] *= utp / u[gtp]

    # determine size of the derivative discontinuity at gtp
    dudrminus = 0.5 * (u[gtp + 1] - u[gtp - 1]) / dr[gtp]
    A = (dudrplus - dudrminus) * utp

    return nodes, A

def shoot_confined(u, l, vr, e, r2dvdr, r, dr, c10, c2, scalarrel=False,
                   gmax=None, rc=10., beta=7.):
    """This method is used by the solve_confined method."""
    # XXX much of this is pasted from the ordinary shoot method

    if scalarrel:
        x = 0.5 * alpha**2 # x = 1 / (2c^2)
        Mr = r * (1.0 + x * e) - x * vr
    else:
        Mr = r
    c0 = l * (l + 1) + 2 * Mr * (vr - e * r)
    if gmax is None and npy.alltrue(c0 > 0):
        print """
Problem with initial electron density guess!  Try to run the program
with the '-n' option (non-scalar-relativistic calculation) and then
try again without the '-n' option (this will generate a good initial
guess for the density).
"""
        raise SystemExit
    c1 = c10
    if scalarrel:
        c0 += x * r2dvdr / Mr
        c1 = c10 - x * r * r2dvdr / (Mr * dr)

    # vectors needed for numeric integration of diff. equation
    fm = 0.5 * c1 - c2
    fp = 0.5 * c1 + c2
    f0 = c0 - 2 * c2

    if gmax is None:
        gcut = int(rc * len(r) / (beta + rc))
        # set boundary conditions at r -> oo (u(oo) = 0 is implicit)
        u[gcut-1] = 1.
        u[gcut:] = 0.

        # perform backwards integration from infinity to the turning point
        g = gcut-2
        u[g] = u[g+1] * f0[g+1] / fm[g+1]
        
        while c0[g] > 0.0: # this defines the classical turning point
            u[g - 1] = (f0[g] * u[g] + fp[g] * u[g + 1]) / fm[g]
            if u[g - 1] < 0.0:
                # There should't be a node here!  Use a more negative
                # eigenvalue:
                print '!!!!!!',
                return 100, None
            if u[g - 1] > 1e100:
                u *= 1e-100
            g -= 1

        # stored values of the wavefunction and the first derivative
        # at the turning point
        gtp = g + 1
        utp = u[gtp]
        dudrplus = 0.5 * (u[gtp + 1] - u[gtp - 1]) / dr[gtp]
    else:
        gtp = gmax

    # set boundary conditions at r -> 0
    u[0] = 0.0
    u[1] = 1.0

    # perform forward integration from zero to the turning point
    g = 1
    nodes = 0
    while g <= gtp: # integrate one step further than gtp
                    # (such that dudr is defined in gtp)
        u[g + 1] = (fm[g] * u[g - 1] - f0[g] * u[g]) / fp[g]
        if u[g + 1] * u[g] < 0:
            nodes += 1
        g += 1
    if gmax is not None:
        return

    # scale first part of wavefunction, such that it is continuous at gtp
    u[:gtp + 2] *= utp / u[gtp]

    # determine size of the derivative discontinuity at gtp
    dudrminus = 0.5 * (u[gtp + 1] - u[gtp - 1]) / dr[gtp]
    A = (dudrplus - dudrminus) * utp

    return nodes, A

if __name__ == '__main__':
    a = AllElectron('Cu', scalarrel=True)
    a.run()

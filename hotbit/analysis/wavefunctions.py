import box
import numpy as nu
import math
import os, pickle
from box import mix
from time import time
from ase import *
acos=math.acos
cos=math.cos
sin=math.sin
atan=math.atan
sqrt=math.sqrt
pi=math.pi
exp=nu.exp



def phival(x,y):
    """ Return azimuthal angle for ALL x,y. """
    e=1E-16
    if x>e and y>e:
        return atan(y/x)
    elif x<-e and y>e:
        return atan(y/x) + pi
    elif x<-e and y<-e:
        return atan(y/x) + pi
    elif x>e and y<-e:
        return atan(y/x)+2*pi
    elif abs(x)<=e and abs(y)<=e:
        return 0.0
    elif x>e and abs(y)<=e:
        return 0.0
    elif y>e and abs(x)<=e:
        return pi/2
    elif x<-e and abs(y)<=e:
        return pi
    elif y<-e and abs(x)<=e:
        return 3*pi/2
    else:
        raise RuntimeError('Strange things in phival')


def to_spherical_coordinates(vec):
    """ Transforms the given cartesien vector to spherical
        coordinates. """
    x, y, z = vec
    r = float(nu.linalg.norm(vec))
    if r < 1e-6:
        r = 0.0
        theta = 0.0
        phi = 0.0
    else:
        theta = float(acos(z/r))
        phi = float(phival(x,y))
    assert 0 <= theta <= nu.pi
    assert 0 <= phi <= 2*nu.pi
    return (r, theta, phi)


def factorial(n):
    if n == 0:
        return 1
    else:
        return n*factorial(n-1)


def Ylm(l, m):
    """ Return the spherical harmonic function. """
    if l == 0:
        if m == 0:
            def ret((theta, phi)):
                return 0.5*sqrt(1/nu.pi)
    if l == 1:
        if m == -1:
            def ret((theta, phi)):
                return 0.5*sqrt(3./(2*pi))*exp(-1j*phi)*sin(theta)
        if m == 0:
            def ret((theta, phi)):
                return 0.5*sqrt(3/pi)*cos(theta)
        if m == 1:
            def ret((theta, phi)):
                return -0.5*sqrt(3./(2*pi))*exp(1j*phi)*sin(theta)
    if l == 2:
        if m == -2:
            def ret((theta, phi)):
                return 0.25*sqrt(15/(2*pi))*exp(-2j*phi)*sin(theta)**2
        if m ==  -1:
            def ret((theta, phi)):
                return 0.5*sqrt(15/(2*pi))*exp(-1j*phi)*sin(theta)*cos(theta)
        if m == 0:
            def ret((theta, phi)):
                return 0.25*sqrt(5/pi)*(3*cos(theta)**2 - 1)
        if m == 1:
            def ret((theta, phi)):
                return -0.5*sqrt(15/(2*pi))*exp(1j*phi)*sin(theta)*cos(theta)
        if m == 2:
            def ret((theta, phi)):
                return 0.25*sqrt(15/(2*pi))*exp(2j*phi)*sin(theta)**2
    if l == 3:
        if m == -3:
            def ret((theta, phi)):
                return 0.125*sqrt(35/pi)*exp(-3j*phi)*sin(theta)**3
        if m == -2:
            def ret((theta, phi)):
                return 0.25*sqrt(105/(2*pi))*exp(-2j*phi)*sin(theta)**2 * cos(theta)
        if m == -1:
            def ret((theta, phi)):
                return 0.125*sqrt(21/pi)*exp(-1j*phi)*sin(theta)*(5*cos(theta)**2 - 1)
        if m == 0:
            def ret((theta, phi)):
                return 0.25*sqrt(7/pi)*(5*cos(theta)**3 - 3*cos(theta))
        if m == 1:
            def ret((theta, phi)):
                return -0.125*sqrt(21/pi)*exp(1j*phi)*sin(theta)*(5*cos(theta)**2 - 1)
        if m == 2:
            def ret((theta, phi)):
                return 0.25*sqrt(105/(2*pi))*exp(2j*phi)*sin(theta)**2 * cos(theta)
        if m == 3:
            def ret((theta, phi)):
                return -0.125*sqrt(35/pi)*exp(3j*phi)*sin(theta)**3
    if l == 4:
        if abs(m) == 4:
            def ret((theta, phi)):
                return (3./16)*sqrt(35/(2*pi))*exp(m*1j*phi)*sin(theta)**4
        if abs(m) == 3:
            def ret((theta, phi)):
                return nu.sign(-m)*(3./8)*sqrt(35/pi)*exp(m*1j*phi)*sin(theta)**3 * cos(theta)
        if abs(m) == 2:
            def ret((theta, phi)):
                return (3./8)*sqrt(5/(2*pi))*exp(m*1j*phi)*sin(theta)**2 * (7*cos(theta)**2 - 1)
        if abs(m) == 1:
            def ret((theta, phi)):
                return nu.sign(-m)*(3./8)*sqrt(5/pi)*exp(m*1j*phi)*sin(theta)*(7*cos(theta)**3 - 3*cos(theta))
        if m == 0:
            def ret((theta, phi)):
                return (3./16)*sqrt(1/pi)*(35*cos(theta)**4 - 30*cos(theta)**2 + 3)

    if l == 5:
        if abs(m) == 5:
            def ret((theta, phi)):
                return nu.sign(-m)*(3./32)*sqrt(77/pi)*exp(m*1j*phi)*sin(theta)**5
        if abs(m) == 4:
            def ret((theta, phi)):
                return (3./16)*sqrt(385/(2*pi))*exp(m*1j*phi)*sin(theta)**4 * cos(theta)
        if abs(m) == 3:
            def ret((theta, phi)):
                return nu.sign(-m)*(1./32)*sqrt(385/pi)*exp(m*1j*phi)*sin(theta)**3 * (9*cos(theta)**2 - 1)
        if abs(m) == 2:
            def ret((theta, phi)):
                return 0.125*sqrt(1155/(2*pi))*exp(m*1j*phi)*sin(theta)**2 * (3*cos(theta)**3 - cos(theta))
        if abs(m) == 1:
            def ret((theta, phi)):
                return nu.sign(-m)*(1./16)*sqrt(165/(2*pi))*exp(m*1j*phi)*sin(theta) * (21*cos(theta)**4 - 14*cos(theta)**2 + 1)
        if m == 0:
            def ret((theta, phi)):
                return (1./16)*sqrt(11/pi) * (63*cos(theta)**5 - 70*cos(theta)**3 + 15*cos(theta))
    if l == 6:
        if abs(m) == 6:
            def ret((theta, phi)):
                return (1./64)*sqrt(3003/pi)*exp(m*1j*phi)*sin(theta)**6
        if abs(m) == 5:
            def ret((theta, phi)):
                return nu.sign(-m)*(3./32)*sqrt(1001/pi)*exp(m*1j*phi)*sin(theta)**5 * cos(theta)
        if abs(m) == 4:
            def ret((theta, phi)):
                return (3./32)*sqrt(91/(2*pi))*exp(m*1j*phi)*sin(theta)**4 * (11*cos(theta)**2 - 1)
        if abs(m) == 3:
            def ret((theta, phi)):
                return nu.sign(-m)*(1./32)*sqrt(1365/pi)*exp(m*1j*phi)*sin(theta)**3 * (11*cos(theta)**3 - 3*cos(theta))
        if abs(m) == 2:
            def ret((theta, phi)):
                return (1./64)*sqrt(1365/pi)*exp(m*1j*phi)*sin(theta)**2 * (33*cos(theta)**4 - 18*cos(theta)**2 + 1)
        if abs(m) == 1:
            def ret((theta, phi)):
                return nu.sign(-m)*(1./16)*sqrt(273/(2*pi))*exp(m*1j*phi)*sin(theta) * (33*cos(theta)**5 - 30*cos(theta)**3 + 5*cos(theta))
        if m == 0:
            def ret((theta, phi)):
                return (1./32)*sqrt(13/pi)*(231*cos(theta)**6 - 315*cos(theta)**4 + 105*cos(theta)**2 - 5)
    return ret


states=['s','px','py','pz','dxy','dyz','dzx','dx2-y2','d3z2-r2']
def angular(r,wf):
    """ Return angular part of wave function.
    
    The real combinations of Y_lm's
    
    parameters:
    -----------
    r: (not normalized) position vector
    wf: index or symbol for state (look below)
    """
    R=sqrt(sum(r**2))
    if R<1E-14:
        return 0.0
    theta=acos(r[2]/R)
    phi=phival(r[0],r[1])

    if type(wf)!=type(1):
        wf=states.index(wf)

    if wf==0:   return 1/sqrt(4*pi)
    elif wf==1: return sqrt(3/(4*pi))*sin(theta)*cos(phi)
    elif wf==2: return sqrt(3/(4*pi))*sin(theta)*sin(phi)
    elif wf==3: return sqrt(3/(4*pi))*cos(theta)
    elif wf==4: return sqrt(15/(4*pi))*sin(theta)**2*cos(phi)*sin(phi)
    elif wf==5: return sqrt(15/(4*pi))*sin(theta)*cos(theta)*sin(phi)
    elif wf==6: return sqrt(15/(4*pi))*sin(theta)*cos(theta)*cos(phi)
    elif wf==7: return 0.5*sqrt(15/(4*pi))*sin(theta)**2*cos(2*phi)
    elif wf==8: return 0.5*sqrt(5/(4*pi))*(3*cos(theta)**2-1)


class WaveFunctions:
    def __init__(self,atoms,dr=0.5):
        """ Wave functions into real grid.
        
        parameters:
        -----------
        dr: (approximate) grid spacing in Angstroms.
        """
        calc=atoms.get_calculator()
        self.atoms=atoms
        self.calc=calc
        self.el=calc.el
        self.st=calc.st
        self.L=self.atoms.get_cell().diagonal() / Bohr
        self.N=ceil(self.L / (dr/Bohr) )
        self.grid=[]
        for i in range(3):
            self.grid.append( nu.linspace(0,self.L[i],self.N[i]) )
        self.dr=self.L/self.N
        self.dV=prod(self.dr)


    def get_wf(self,i):
        """ Return wave function i on given grid. """
        wf=self.st.wf[:,i].copy()
        orbs=self.el.orbitals()
        wfg=nu.zeros(self.N,type(wf))
        for orb,c in zip(orbs,wf):
            symb, orbtype, Rnl=orb['symbol'], orb['orbital'], orb['Rnl']
            assert abs(Rnl(5.0))>1E-9
            for i,x in enumerate(self.grid[0]):
                for j,y in enumerate(self.grid[1]):
                    for k,z in enumerate(self.grid[2]):
                        r0=nu.array([x,y,z])
                        r=self.el.vector(orb['atom'],rj=r0)
                        wfg[i,j,k] += c*Rnl(mix.norm(r))*angular(r,orbtype)
        return wfg


    def get_wf_pseudo_density(self,i):
        """ Return the density from localized pseudo wave functions for state i."""
        rho=abs(self.get_wf(i))**2
        rho/=float(rho.flatten().sum()*self.dV) #ensure normalization
        return rho


    def get_pseudo_density(self):
        """ Return the valence electron density from pseudo wave functions."""
        rho=None
        for i,f in enumerate(self.st.f):
            if f<1E-4: break
            rhoi = self.get_wf_pseudo_density(i)
            if rho==None: rho = rhoi*float(f)
            else:
                rho+=rhoi*float(f)
        return rho


        rho=self.get_wf(i)
        rho*=rho.conjugate()
        rho/=rho.flatten().sum() #ensure normalization
        return rho


    def write_vtk(self,i,fname=None):
        """ Write .vtk file of wave function with *index* i. """
        wf=self.st.wf[:,i].copy()
        orbs=self.el.orbitals()
        wfg=nu.zeros(self.N)
        grid=[]
        for i in range(3):
            grid.append( nu.linspace(0,self.L[i],self.N[i]) )
        for orb,c in zip(orbs,wf):
            symb, orbtype, Rnl=orb['symbol'], orb['orbital'], orb['Rnl']
            for i,x in enumerate(grid[0]):
                for j,y in enumerate(grid[1]):
                    for k,z in enumerate(grid[2]):
                        r0=nu.array([x,y,z])
                        r=self.el.vector(orb['atom'],rj=r0)
                        wfg[i,j,k]+=c*Rnl(mix.norm(r))*angular(r,orbtype)
        box.vtk.rectilinear_vtk(grid,wfg,fname)


class JelliumAnalysis:


    def __init__(self, atoms, origin=None, maxl=3, R_0=3, a=0.2):
        """
        A class to analyse the wave functions of the system by projecting
        them onto the spherical harmonics.

        atoms:  The ase atoms object with a ground state calculated
                of the expansion data file calculated before.
        origin: The center of the expansion (Ang). If None, the center of
                the mass will be used.
        maxl:   The largest angular momentum the expansion is performed.
        R_0:    The radius of the expansion (Ang)
        a:      The approximate length of the side of cubic grid box (Ang)
        N:      Number of grid points in expansion sphere in one dimension
                (the total number is then N^3).
        """
        # values that are saved/loaded
        self.data = ['R_0',
                     'origin',
                     'N',
                     'a',
                     'dV',
                     'grid_points',
                     'maxl',
                     'l_array',
                     'e',
                     'occ',
                     'fermi_level',
                     'c_nl',
                     'weights',
                     'dim',
                     'N_wf',
                     'needed_orbitals']
        self.letters = "spdfghi"
        self.loaded = None

        if type(atoms) == str:
            self.load(atoms)
        else:
            self.calc = atoms.get_calculator()
            self.R_0 = R_0 / Bohr

            self.atoms = atoms
            if origin == None:
                self.origin = self.calc.el.get_center_of_mass()
            else:
                self.origin = nu.array(origin)/Bohr
            self.maxl = maxl
            self.l_array = range(min(7, maxl+1))
            self.N_wf = self.calc.st.norb
            self.fermi_level = self.calc.get_fermi_level()
            self.e = self.calc.st.get_eigenvalues() * Hartree
            self.occ = self.calc.st.get_occupations()

            # the expansion coefficients
            self.c_nl = nu.zeros((self.calc.st.norb, len(self.l_array)))

            # The amount of states inside the expansion radius
            self.weights = nu.zeros(self.calc.st.norb)
            self.basis_functions = {}

            self.create_uniform_cubic_grid(a)
            self.mark_grids()
            self.mark_needed_basis_functions()
        self.greetings()


    def load(self, file):
        """ Load expansion data to file. """
        f = open(file)
        while True:
            try:
                name, data = pickle.load(f)
                self.__dict__[name] = data
            except EOFError:
                break
        self.loaded = file
        f.close()


    def write(self, file):
        """ Write expansion data to file. """
        f = open(file, 'w')
        for name in self.data:
            pickle.dump([name, self.__dict__[name]], f)
        f.close()


    def create_uniform_cubic_grid(self, a):
        """
        Creates a grid of N x N x N cubes. The grid is built so that
        the origin of the expansion is between the grid points.
        """
        a = a / Bohr # the requested grid spacing
        self.N = nu.ceil(2*self.R_0 / a)
        # we want an odd number of grid points per side, so that the
        # center of the expansion is one grid point.
        if self.N % 2 == 0:
            self.N += 1
        self.a = 2*self.R_0/self.N # the actual grid spacing
        assert abs(self.N*self.a - 2*self.R_0) < 1e-3
        self.dV = self.a**3
        self.grid_points = []
        for d in range(3):
            R = self.origin[d]
            g = nu.linspace(R - self.R_0 + self.a/2, R + self.R_0 - self.a/2, self.N)
            self.grid_points.append(g)
        self.dim = nu.array([len(axis) for axis in self.grid_points])


    def estimate_memory_consumption(self, human_readable=False):
        """ Give an approximation on how much the arrays need memory. """
        # grid points for one state
        N = nu.prod(self.dim)
        #n_orb = self.calc.el.norb
        n_orb = nu.sum(nu.where(self.needed_orbitals, 1, 0))
        mem = 0.0
        # the spherical harmonics
        mem += 0.5 * self.maxl * (self.maxl+1) * nu.array([1], nu.complex).itemsize * N
        # the basis functions
        mem += n_orb * nu.array([1], nu.float).itemsize * N
        # working arrays
        mem += 4 * nu.array([1], nu.float).itemsize * N
        f = 'byte'
        if mem / 1024.0 > 1.0:
            memm = mem / 1024.0
            f = 'kb'
        if memm / 1024.0 > 1.0:
            memm = memm / 1024.0
            f = 'Mb'
        if memm / 1024.0 > 1.0:
            memm = memm / 1024.0
            f = 'Gb'
        if human_readable:
            return "%0.3f %s" % (memm, f)
        else:
            return mem


    def mark_grids(self):
        """ Create a grid that contains indices that tell to which
            shell the grid belongs to.
              Index = 0: outside the expansion shell
              Index >= 1: the index of the shell (1 is the most inner one)
            The thickness of the shells is a. """
        self.shell_index_grid = nu.zeros(self.dim, dtype=int)
        self.shells = {}
        for i, x in enumerate(self.grid_points[0]):
            for j, y in enumerate(self.grid_points[1]):
                for k, z in enumerate(self.grid_points[2]):
                    vec = nu.array((x,y,z)) - self.origin
                    norm = nu.linalg.norm(vec)
                    if norm <= self.R_0:
                        index = int(round(norm/self.a)) + 1
                        if index in self.shells:
                            self.shells[index] += 1
                        else:
                            self.shells[index] = 1
                        self.shell_index_grid[i,j,k] = index


    def mark_needed_basis_functions(self):
        """ Chech which basis functions are needed in order to calculate
            the wave functions inside the expansion sphere. """
        self.needed_orbitals = nu.zeros(self.calc.st.norb ,dtype=bool)
        positions = self.calc.el.get_positions()
        orbitals = self.calc.el.orbitals()
        for m in range(self.calc.st.norb):
            orbital = orbitals[m]
            atom = orbital['atom']
            orbital_type = orbital['orbital']
            symbol = orbital['symbol']
            wf_range = self.calc.el.elements[symbol].get_wf_range(orbital_type)
            if wf_range == None:
                raise Exception('No radial function %s found for %s (maybe the element file does not contain it?).' % (orbital_type, symbol))
            R = positions[atom]
            if nu.linalg.norm(R - self.origin) < self.R_0 + wf_range:
                self.needed_orbitals[m] = True


    def ylms_to_grid(self):
        """ Calculates the values of spherical harmonics centered
            to the origin of the expansion to the grid. """
        self.ylms = {}
        origin = self.origin
        for l in self.l_array:
            for m in range(-l,l+1):
                y_lm = Ylm(l, m)
                t1 = time()
                print "Calculating Y_(l=%i, m=%s) to grid..." % (l, repr(m).rjust(2)),
                values = nu.zeros(self.dim, dtype=nu.complex)
                for i, x in enumerate(self.grid_points[0]):
                    for j, y in enumerate(self.grid_points[1]):
                        for k, z in enumerate(self.grid_points[2]):
                            r, theta, phi = to_spherical_coordinates(nu.array((x,y,z))-self.origin)
                            values[i,j,k] = y_lm((theta, phi))
                print "in %i seconds." % (time() - t1)
                self.ylms[l,m] = values


    def get_ylm(self, l, m):
        """ Return spherical harmonic Y_lm on grid. """
        return self.ylms[l,m]


    def get_basis_function(self, m):
        """ Return m:th basis function on grid. """
        if not m in self.basis_functions:
            orbital = self.calc.el.orbitals()[m]
            atom = orbital['atom']
            R = self.calc.el.get_positions()[atom]
            orbital_type = orbital['orbital']
            if self.needed_orbitals[m] == False:
                print "  Basis function %i/%i is outside the range." % (orbital['index']+1, self.calc.st.norb)
                self.basis_functions[m] = None
            else:
                t1 = time()
                print "  Calculating basis function %i/%i to grid..." % (orbital['index']+1, self.calc.st.norb),
                r_nl = orbital['Rnl']
                basis_function = nu.zeros(self.dim, dtype=nu.float)
                for i, x in enumerate(self.grid_points[0]):
                    for j, y in enumerate(self.grid_points[1]):
                        for k, z in enumerate(self.grid_points[2]):
                            vec = nu.array((x,y,z))-nu.array(R)
                            r, theta, phi = to_spherical_coordinates(vec)
                            basis_function[i,j,k] = r_nl(r)*angular(vec, orbital_type)
                print "  in %i seconds." % (time() - t1)
                self.basis_functions[m] = basis_function
        return self.basis_functions[m]


    def get_state(self, k):
        """ Return the k:th wave function inside the expansion grid. """
        wf_coeffs = self.calc.st.wf[:,k]
        state_grid = nu.zeros(self.dim, dtype=nu.float)
        t1 = time()
        for wf_coef, orb in zip(wf_coeffs, self.calc.el.orbitals()):
            basis_function = self.get_basis_function(orb['index'])
            if basis_function != None:
                state_grid += wf_coef * basis_function
        return state_grid


    def analyse_states(self):
        """ Perform the angular momentum analysis on all states. """
        for n in range(self.calc.st.norb):
            print "Analysing state %i..." % n
            t1 = time()
            state_grid = self.get_state(n)
            state_grid_squared = state_grid.conjugate() * state_grid
            self.weights[n] = nu.sum(state_grid_squared * nu.where(self.shell_index_grid != 0, 1, 0)) * self.dV
            for l in self.l_array:
                self.c_nl[n,l] = self.weights[n] * self.analyse_state(state_grid, l)
            print "State %i analyzed in %i seconds." % (n, time() -  t1)
        for n in range(self.calc.st.norb):
            if float(nu.sum(self.c_nl[n,:])) > 0.01:
                self.c_nl[n,:] = self.c_nl[n,:]/float(nu.sum(self.c_nl[n,:]))


    def analyse_state(self, state_grid, l):
        """ Performs the angular momentum analysis with respect to
            angular momentum l to the state n. """
        c = 0.0
        for m in range(-l,l+1):
            ylm = self.get_ylm(l, m)
            # The integration
            for i in self.shells.keys():
                # the mask that gives the grid points of the shell
                shell_grid = nu.where(self.shell_index_grid == i, 1, 0)
                # the number of boxes in the i:th shell
                N_shell = self.shells[i]
                # the integration over the solid angle
                phi_nlm = nu.sum(shell_grid * ylm.conjugate() * state_grid)
                c += phi_nlm.conjugate() * phi_nlm * 4*pi*self.dV / N_shell
        return c


    def greetings(self):
        print "\n*** Starting the angular momentum analysis. ***"
        if self.loaded != None:
            print "Using data from %s" % self.loaded
        print "The grid contains %i x %i x %i grid points" % tuple(self.dim)
        print "The grid spacing is %0.3f (Ang)" % (self.a*Bohr)
        print "The center of the expansion (Ang): %0.2f, %0.2f, %0.2f" % tuple(self.origin * Bohr)
        print "The radius of the expansion: %0.2f (ang)" % (self.R_0 * Bohr)
        print "The analysis is performed on angular momenta:",
        for l in self.l_array:
            print self.letters[l],
        print ""
        if self.loaded == None:
            print "Estimated amount of memory required: %s" % self.estimate_memory_consumption(human_readable=True)
        print ""


    def write_readable(self, file):
        """ Write the spherical harmonic expansion coefficients to file. """
        f = open(file, 'w')
        print >> f, "# The center of the expansion (in Ang): %0.2f, %0.2f, %0.2f" % tuple(self.origin * Bohr)
        print >> f, "# The radius of the expansion (in Ang): %0.2f" % (self.R_0 * Bohr)
        print >> f, "# The shell thickness (in Ang): %0.2f" % (self.a * Bohr)
        print >> f, "#state  energy(eV)   weight     occ",
        for l in self.l_array:
            print >> f, "%6s" % self.letters[l],
        print >> f, ""
        for n in range(self.calc.st.norb):
            print >> f, "%5i %12.4f %7.4f %7.4f" % (n, e[n], w[n], occ[n]),
            for l in self.l_array:
                print >> f, "%6.3f" % self.c_nl[n,l],
            print >> f, ""
        f.close()


    def gaussian_peak(self, x, x0, width):
        return exp( - (x-x0)**2 / (4*width**2) )


    def make_plot(self, width=0.1, bands=None, limits=None, out=None):
        import pylab
        colors = ['#FFFF00','#FF0000','#2758D3','#5FD300',
                  '#058C00','#E1AB18','#50E1D0']

        if bands == None:
            bands = self.N_wf
        c_nl = self.c_nl
        e = self.e[:bands] - self.fermi_level
        if limits == None:
            e_min, e_max = min(e), max(e)
            empty = 0.1*(e_max - e_min)
        else:
            e_min, e_max = limits
            empty = 0.0
        x_grid = nu.linspace(e_min-empty, e_max+empty, 2000)
        y_0 = nu.zeros(len(x_grid))
        for l in self.l_array:
            y_1 = y_0.copy()
            for n, eig in enumerate(e):
                y_1 += c_nl[n,l]*self.gaussian_peak(x_grid, eig, width)
            x, y = pylab.poly_between(x_grid, y_0, y_1)
            pylab.fill(x, y, facecolor=colors[l], edgecolor='none', label=self.letters[l])
            y_0 = y_1
        pylab.xlabel("Energy (eV)")
        pylab.ylabel("Arbitrary units")
        pylab.title("Angular momentum analysis")
        pylab.legend()
        if out == None:
            out = 'Jellium_analysis_w%0.3f.eps'
        pylab.savefig(out)
        pylab.show()


    def run(self):
        if self.loaded == None:
            self.ylms_to_grid()
            self.analyse_states()
            self.write("ylm_expansion.hb")



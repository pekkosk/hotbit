import box
import numpy as np
import math
import os, pickle
from box import mix
import time
from ase import *
from box.broaden import *
from box.mix import phival
acos=math.acos
cos=math.cos
sin=math.sin
atan=math.atan
sqrt=math.sqrt
pi=math.pi
exp=np.exp




def to_spherical_coordinates(vec):
    """ Transforms the given cartesien vector to spherical
        coordinates. """
    x, y, z = vec
    r = float(np.linalg.norm(vec))
    if r < 1e-6:
        r = 0.0
        theta = 0.0
        phi = 0.0
    else:
        theta = float(acos(z/r))
        phi = float(phival(x,y))
    assert 0 <= theta <= np.pi
    assert 0 <= phi <= 2*np.pi
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
                return 0.5*sqrt(1/np.pi)
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
                return np.sign(-m)*(3./8)*sqrt(35/pi)*exp(m*1j*phi)*sin(theta)**3 * cos(theta)
        if abs(m) == 2:
            def ret((theta, phi)):
                return (3./8)*sqrt(5/(2*pi))*exp(m*1j*phi)*sin(theta)**2 * (7*cos(theta)**2 - 1)
        if abs(m) == 1:
            def ret((theta, phi)):
                return np.sign(-m)*(3./8)*sqrt(5/pi)*exp(m*1j*phi)*sin(theta)*(7*cos(theta)**3 - 3*cos(theta))
        if m == 0:
            def ret((theta, phi)):
                return (3./16)*sqrt(1/pi)*(35*cos(theta)**4 - 30*cos(theta)**2 + 3)

    if l == 5:
        if abs(m) == 5:
            def ret((theta, phi)):
                return np.sign(-m)*(3./32)*sqrt(77/pi)*exp(m*1j*phi)*sin(theta)**5
        if abs(m) == 4:
            def ret((theta, phi)):
                return (3./16)*sqrt(385/(2*pi))*exp(m*1j*phi)*sin(theta)**4 * cos(theta)
        if abs(m) == 3:
            def ret((theta, phi)):
                return np.sign(-m)*(1./32)*sqrt(385/pi)*exp(m*1j*phi)*sin(theta)**3 * (9*cos(theta)**2 - 1)
        if abs(m) == 2:
            def ret((theta, phi)):
                return 0.125*sqrt(1155/(2*pi))*exp(m*1j*phi)*sin(theta)**2 * (3*cos(theta)**3 - cos(theta))
        if abs(m) == 1:
            def ret((theta, phi)):
                return np.sign(-m)*(1./16)*sqrt(165/(2*pi))*exp(m*1j*phi)*sin(theta) * (21*cos(theta)**4 - 14*cos(theta)**2 + 1)
        if m == 0:
            def ret((theta, phi)):
                return (1./16)*sqrt(11/pi) * (63*cos(theta)**5 - 70*cos(theta)**3 + 15*cos(theta))
    if l == 6:
        if abs(m) == 6:
            def ret((theta, phi)):
                return (1./64)*sqrt(3003/pi)*exp(m*1j*phi)*sin(theta)**6
        if abs(m) == 5:
            def ret((theta, phi)):
                return np.sign(-m)*(3./32)*sqrt(1001/pi)*exp(m*1j*phi)*sin(theta)**5 * cos(theta)
        if abs(m) == 4:
            def ret((theta, phi)):
                return (3./32)*sqrt(91/(2*pi))*exp(m*1j*phi)*sin(theta)**4 * (11*cos(theta)**2 - 1)
        if abs(m) == 3:
            def ret((theta, phi)):
                return np.sign(-m)*(1./32)*sqrt(1365/pi)*exp(m*1j*phi)*sin(theta)**3 * (11*cos(theta)**3 - 3*cos(theta))
        if abs(m) == 2:
            def ret((theta, phi)):
                return (1./64)*sqrt(1365/pi)*exp(m*1j*phi)*sin(theta)**2 * (33*cos(theta)**4 - 18*cos(theta)**2 + 1)
        if abs(m) == 1:
            def ret((theta, phi)):
                return np.sign(-m)*(1./16)*sqrt(273/(2*pi))*exp(m*1j*phi)*sin(theta) * (33*cos(theta)**5 - 30*cos(theta)**3 + 5*cos(theta))
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




    # def write_vtk(self,i,fname=None): ----------------------------------------
        # """ Write .vtk file of wave function with *index* i. """ -------------
        # wf=self.st.wf[i,:].copy() --------------------------------------------
        # orbs=self.el.orbitals() ----------------------------------------------
        # wfg=np.zeros(self.N) -------------------------------------------------
        # grid=[] --------------------------------------------------------------
        # for i in range(3): ---------------------------------------------------
            # grid.append( np.linspace(0,self.L[i],self.N[i]) ) ----------------
        # for orb,c in zip(orbs,wf): -------------------------------------------
            # symb, orbtype, Rnl=orb['symbol'], orb['orbital'], orb['Rnl'] -----
            # for i,x in enumerate(grid[0]): -----------------------------------
                # for j,y in enumerate(grid[1]): -------------------------------
                    # for k,z in enumerate(grid[2]): ---------------------------
                        # r0=np.array([x,y,z]) ---------------------------------
                        # r=self.el.vector(orb['atom'],rj=r0) ------------------
                        # wfg[i,j,k]+=c*Rnl(mix.norm(r))*angular(r,orbtype) ----
        # box.vtk.rectilinear_vtk(grid,wfg,fname) ------------------------------


class JelliumAnalysis:


    def __init__(self, atoms, origin=None, maxl=3, R_0=3, a=0.2, file='ylm_expansion.hb'):
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
        file:   The filename where the ylm-expansion data is saved.
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
                     'finished',
                     'c_nl',
                     'weights',
                     'dim',
                     'norb',
                     'analyzed_states',
                     'needed_orbitals',
                     'norb_needed',
                     'file',
                     'letters']
        self.letters = "spdfghi"[0:maxl+1]
        self.loaded = False
        if os.path.isfile(file):
            self.load(file)
        else:
            atoms.get_potential_energy()
            calc = atoms.get_calculator()
            self.R_0 = R_0 / Bohr
            self.finished = False
            if origin == None:
                self.origin = calc.el.get_center_of_mass()
            else:
                self.origin = np.array(origin)/Bohr
            self.create_uniform_cubic_grid(a)
            self.maxl = maxl
            self.l_array = range(min(7, maxl+1))
            self.norb = calc.st.norb
            self.fermi_level = calc.st.occu.get_mu()
            self.e = calc.st.get_eigenvalues()
            self.occ = calc.st.get_occupations()
            self.file = file
            self.c_nl = np.zeros((self.norb, len(self.l_array)))
            self.weights = np.zeros(self.norb)
            self.analyzed_states = np.zeros(self.norb, dtype=bool)

        self.atoms = atoms
        self.calc = atoms.get_calculator()

        self.basis_functions = {}

        self.mark_grids()
        self.mark_needed_basis_functions()

        self.log = open('jellium_analysis.log','a')
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
        self.loaded = True
        f.close()


    def write(self):
        """ Write expansion data to file. """
        f = open(self.file, 'w')
        for name in self.data:
            pickle.dump([name, self.__dict__[name]], f)
        f.close()


    def create_uniform_cubic_grid(self, a):
        """
        Creates a grid of N x N x N cubes. The grid is built so that
        the origin of the expansion is between the grid points.
        """
        a = a / Bohr # the requested grid spacing
        self.N = np.ceil(2*self.R_0 / a)
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
            g = np.linspace(R - self.R_0 + self.a/2, R + self.R_0 - self.a/2, self.N)
            self.grid_points.append(g)
        self.dim = np.array([len(axis) for axis in self.grid_points])


    def estimate_memory_consumption(self, human_readable=False):
        """ Give an approximation on how much the arrays need memory. """
        # grid points for one state
        N = np.prod(self.dim)
        mem = 0.0
        # the spherical harmonics
        mem += (self.maxl+1)**2 * np.array([1], np.complex).itemsize * N
        # the basis functions
        mem += self.norb_needed * np.array([1], np.float).itemsize * N
        # working arrays
        mem += 4 * np.array([1], np.float).itemsize * N
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
        self.shell_index_grid = np.zeros(self.dim, dtype=int)
        self.shells = {}
        for i, x in enumerate(self.grid_points[0]):
            for j, y in enumerate(self.grid_points[1]):
                for k, z in enumerate(self.grid_points[2]):
                    vec = np.array((x,y,z)) - self.origin
                    norm = np.linalg.norm(vec)
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
        self.needed_orbitals = np.zeros(self.norb ,dtype=bool)
        positions = self.calc.el.get_positions()
        orbitals = self.calc.el.orbitals()
        for m in range(self.norb):
            orbital = orbitals[m]
            atom = orbital['atom']
            orbital_type = orbital['orbital']
            symbol = orbital['symbol']
            wf_range = self.calc.el.elements[symbol].get_wf_range(orbital_type)
            if wf_range == None:
                raise Exception('No radial function %s found for %s (maybe the element file does not contain it?).' % (orbital_type, symbol))
            R = positions[atom]
            if np.linalg.norm(R - self.origin) < self.R_0 + wf_range:
                self.needed_orbitals[m] = True
        self.norb_needed = np.sum(np.where(self.needed_orbitals, 1, 0))


    def ylms_to_grid(self):
        """ Calculates the values of spherical harmonics centered
            to the origin of the expansion to the grid. """
        self.ylms = {}
        origin = self.origin
        for l in self.l_array:
            for m in range(-l,l+1):
                y_lm = Ylm(l, m)
                t1 = time.time()
                print >> self.log, "Calculating Y_(l=%i, m=%s) to grid..." % (l, repr(m).rjust(2))
                self.log.flush()
                values = np.zeros(self.dim, dtype=np.complex)
                for i, x in enumerate(self.grid_points[0]):
                    for j, y in enumerate(self.grid_points[1]):
                        for k, z in enumerate(self.grid_points[2]):
                            r, theta, phi = to_spherical_coordinates(np.array((x,y,z))-self.origin)
                            values[i,j,k] = y_lm((theta, phi))
                print >> self.log, "  in %i seconds." % (time.time() - t1)
                self.ylms[l,m] = values


    def basis_functions_to_grid(self):
        """ Calculate the basis functions into the grid. """
        for m in range(self.norb):
            orbital = self.calc.el.orbitals()[m]
            atom = orbital['atom']
            R = self.calc.el.get_positions()[atom]
            orbital_type = orbital['orbital']
            if self.needed_orbitals[m] == False:
                print >> self.log, "Basis function %i/%i is outside the range." % (orbital['index']+1, self.norb)
                self.basis_functions[m] = None
            else:
                t1 = time.time()
                print >> self.log, "Calculating basis function %i/%i to grid..." % (orbital['index']+1, self.norb)
                self.log.flush()
                r_nl = orbital['Rnl']
                basis_function = np.zeros(self.dim, dtype=np.float)
                for i, x in enumerate(self.grid_points[0]):
                    for j, y in enumerate(self.grid_points[1]):
                        for k, z in enumerate(self.grid_points[2]):
                            vec = np.array((x,y,z))-np.array(R)
                            r, theta, phi = to_spherical_coordinates(vec)
                            basis_function[i,j,k] = r_nl(r)*angular(vec, orbital_type)
                print >> self.log, "  in %i seconds." % (time.time() - t1)
                self.basis_functions[m] = basis_function


    def get_ylm(self, l, m):
        """ Return spherical harmonic Y_lm on grid. """
        return self.ylms[l,m]


    def get_basis_function(self, m):
        """ Return m:th basis function on grid. """
        return self.basis_functions[m]


    def get_state(self, k):
        """ Return the k:th wave function inside the expansion grid. """
        wf_coeffs = self.calc.st.wf[0,k,:]
        assert np.sum(wf_coeffs.imag < 1e-3)
        wf_coeffs = wf_coeffs.real
        state_grid = np.zeros(self.dim, dtype=np.float)
        for wf_coef, orb in zip(wf_coeffs, self.calc.el.orbitals()):
            basis_function = self.get_basis_function(orb['index'])
            if basis_function != None:
                state_grid += wf_coef * basis_function
        return state_grid


    def analyse_states(self):
        """ Perform the angular momentum analysis on all states. """
        for n in range(self.norb):
            if self.analyzed_states[n] == False:
                print >> self.log, "Analysing state no. %i (%i/%i)..." % (n, n+1, self.norb)
                self.log.flush()
                t1 = time.time()
                state_grid = self.get_state(n)
                state_grid_squared = state_grid.conjugate() * state_grid
                self.weights[n] = np.sum(state_grid_squared * np.where(self.shell_index_grid != 0, 1, 0)) * self.dV
                for l in self.l_array:
                    self.c_nl[n,l] = self.weights[n] * self.analyse_state(state_grid, l)
                print >> self.log, "  %i seconds." % (time.time() -  t1)
                self.analyzed_states[n] = True
                self.write()
        self.finished = True
        for n in range(self.norb):
            if float(np.sum(self.c_nl[n,:])) > 0.01:
                self.c_nl[n,:] = self.c_nl[n,:]/float(np.sum(self.c_nl[n,:]))
        self.log.flush()
        self.write_readable('ja.dat')


    def analyse_state(self, state_grid, l):
        """ Performs the angular momentum analysis with respect to
            angular momentum l to the state n. """
        c = 0.0
        for m in range(-l,l+1):
            ylm = self.get_ylm(l, m)
            # The integration
            for i in self.shells.keys():
                # the mask that gives the grid points of the shell
                shell_grid = np.where(self.shell_index_grid == i, 1, 0)
                # the number of boxes in the i:th shell
                N_shell = self.shells[i]
                # the integration over the solid angle
                phi_nlm = np.sum(shell_grid * ylm.conjugate() * state_grid)
                c += phi_nlm.conjugate() * phi_nlm * 4*pi*self.dV / N_shell
        return c


    def greetings(self):
        print >> self.log, "\n*** Starting the angular momentum analysis. ***"
        if self.loaded:
            print >> self.log, "Using data from %s" % self.file
        print >> self.log, "The grid contains %i x %i x %i grid points" % tuple(self.dim)
        print >> self.log, "The grid spacing is %0.3f (Ang)" % (self.a*Bohr)
        print >> self.log, "The center of the expansion (Ang): %0.2f, %0.2f, %0.2f" % tuple(self.origin * Bohr)
        print >> self.log, "The radius of the expansion: %0.2f (ang)" % (self.R_0 * Bohr)
        print >> self.log, "The analysis is performed on angular momenta:",
        for l in self.l_array:
            print >> self.log, self.letters[l],
        print >> self.log, ""
        print >> self.log, "There are %i basis functions, %i are needed." % (self.norb, self.norb_needed)
        print >> self.log, "There are %i/%i states left to analyze." % (np.sum(np.where(self.analyzed_states, 0, 1)), self.norb)
        if self.finished == False:
            print >> self.log, "Estimated amount of memory required: %s" % (self.estimate_memory_consumption(human_readable=True))
        print >> self.log, ""
        self.log.flush()


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
        for n in range(self.norb):
            print >> f, "%5i %12.4f %7.4f %7.4f" % (n, self.e[0,n]*Hartree, self.weights[n], self.occ[0,n]),
            for l in self.l_array:
                print >> f, "%6.3f" % self.c_nl[n,l],
            print >> f, ""
        f.close()


    def run(self):
        if self.finished == False:
            self.atoms.get_potential_energy()
            self.ylms_to_grid()
            self.basis_functions_to_grid()
            self.analyse_states()
            self.write()


def make_plot(file, width=0.1, xlimits=None, ylimits=None, out=None):
    """ Make a cumulative plot from the angular momentum analysis
        of the electron states.
        
        file: data file produced by the Jellium analysis object
        width: the FWHM of the gaussians used to broaden the energies
        xlimits: the limits of the x-axis [xmin, xmax]
        ylimits: the limits of the y-axis [ymin, ymax]
        out: the output file (None=screen)

    """
    import pylab
    f = open(file)
    data = {}
    while True:
        try:
            name, dat = pickle.load(f)
            data[name] = dat
        except EOFError:
            break
    print ""
    print ""
    print "*** Spherical harmonics analysis ***"
    print "center:", data["origin"] * Bohr, "Ang"
    print "radius:", data["R_0"] * Bohr, "Ang"
    print "grid box size (shell thickness):", data["a"] * Bohr, "Ang"
    print "analysis performed on angular momenta:", data["letters"].replace("",",")[1:-1]
    print "Fermi level (subtracted):", data["fermi_level"] * Hartree, "eV"

    pylab.figure()
    colors = ['#FFFF00','#FF0000','#2758D3','#5FD300',
              '#058C00','#E1AB18','#50E1D0']

    weights = data["c_nl"].transpose()
    e = (data["e"][0] - data["fermi_level"]) * Hartree
    if xlimits == None:
        e_min, e_max = min(e), max(e)
        empty = 0.1*(e_max - e_min)
    else:
        e_min, e_max = xlimits
        empty = 0.0
    x_grid = np.linspace(e_min-empty, e_max+empty, 2000)
    make_cumulative_plot(x_grid, e, width, weights,
                         labels=data["letters"], colors=colors)
    pylab.xlabel("Energy (eV)")
    pylab.ylabel("Arbitrary units")
    pylab.title("Angular momentum analysis")
    pylab.legend()
    pylab.ylim(ylimits)
    if out == None:
        pylab.show()
    else:
        pylab.savefig(out)


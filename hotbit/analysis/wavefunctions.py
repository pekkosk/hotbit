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
    def __init__(self,atoms,mesh=(20,20,20)):
        """ Wave functions into real grid object. 
        
        parameters:
        -----------
        mesh: x,y,and z divisions for the grid
        """
        calc=atoms.get_calculator()
        self.atoms=atoms
        self.calc=calc
        self.el=calc.el
        self.st=calc.st
        cell=self.atoms.get_cell()
        self.L=nu.array([cell[0,0],cell[1,1],cell[2,2]])
        self.N=mesh
        
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

    def transform(self, name):
        """ Returns the corresponding angular and magnetic momentum
        indices for the given orbital name. """
        names = ['s','px','py','pz','dxy','dyz','dzx','dx2-y2','d3z2-r2']
        list = ((0,0),(1,-1),(1,0),(1,1),(2,-2),(2,-1),(2,0),(2,1),(2,2))
        return list[names.index(name)]


    def to_spherical_coordinates(self, vec):
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


    def __init__(self, atoms, origin=None, maxl=3, R_0=1, a=0.2, filename='JelliumAnalysis.dat', store_grids=False):
        """
        A class to analyse the wave functions of the system by projecting
        them onto the spherical harmonics.

        origin: The center of the expansion (Ang). If None, the center of
                the mass will be used.
        maxl:   The largest angular momentum the expansion is performed.
        R_0:    The radius of the expansion (Ang)
        a:      The length of the side of cubic grid box (Ang)
        """
        self.a = a/Bohr
        self.calc = atoms.get_calculator()
        self.R_0 = R_0/Bohr
        self.atoms = atoms
        if origin == None:
            self.origin = self.calc.el.get_center_of_mass()
        else:
            self.origin = nu.array(origin)/Bohr
        self.l_array = range(min(7, maxl+1))
        self.filename = filename
        self.store_grids = store_grids
        self.letters = "spdfghi"

        self.c_nl = nu.zeros((self.calc.st.norb, len(self.l_array)))

        # The amount of states inside the expansion radius
        self.weights = nu.zeros(self.calc.st.norb)
        # The norms of the states
        self.norms = nu.zeros(self.calc.st.norb)

        self.create_uniform_cubic_grid()


    def create_uniform_cubic_grid(self):
        """
        Creates a grid of a x a x a cubes, where the grid points
        are in the center of the cubes. The grid is built so that
        the origin of the expansion is one of the grid points.
        """
        a = self.a
        self.dV = a**3
        cell = self.calc.el.get_box_lengths()
        self.grid_points = []
        for d in range(3):
            g = []
            R = self.origin[d]
            limits = [0, cell[d]]
            for i, m in enumerate([-1,1]):
                n = 0
                while True:
                    n += 1
                    if nu.sign(-m)*( R + m*(n-0.5)*a - limits[i] ) > 0:
                        g.append(R + m*n*a)
                    else:
                        break
            g.sort()
            self.grid_points.append(g)
        self.dim = [len(axis) for axis in self.grid_points]


    def mark_grids(self):
        """ Create a grid that contains indices that tell to which
            shell the grid belongs to. The thickness of the shells
            is a. """
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


    def ylms_to_grid(self):
        """ Calculates the values of spherical harmonics centered
            to the origin of the expansion to the grid. """
        self.ylms = {}
        origin = self.origin
        for l in self.l_array:
            for m in range(-l,l+1):
                if os.path.isfile("%s/spherical_harmonic_%i%i.hb" % (os.getcwd(), l,m)):
                    print "found spherical harmonic l=%i,m=%i" % (l, m)
                    f = open("%s/spherical_harmonic_%i%i.hb" % (os.getcwd(), l, m),'r')
                    values = pickle.load(f)
                    f.close()
                    assert nu.values(shape) == self.dim
                else:
                    y_lm = Ylm(l, m)
                    t1 = time()
                    print "Calculating Y_(l=%i,m=%i) to grid..." % (l, m),
                    values = nu.zeros(self.dim, dtype=nu.complex)
                    for i, x in enumerate(self.grid_points[0]):
                        for j, y in enumerate(self.grid_points[1]):
                            for k, z in enumerate(self.grid_points[2]):
                                r, theta, phi = self.to_spherical_coordinates(nu.array((x,y,z))-self.origin)
                                values[i,j,k] = y_lm((theta, phi))
                    print "in %i seconds." % (time() - t1)
                    if self.store_grids:
                        f = open("%s/spherical_harmonic_%i%i.hb" % (os.getcwd(), l, m),'w')
                        pickle.dump(values, f)
                        f.close()
                self.ylms["%i,%i" % (l, m)] = values


    def basis_functions_to_grid(self):
        """ Calculates the basis functions into the cubic grid. """
        self.basis_functions = {}
        orbitals = self.calc.el.orbitals()
        positions = self.calc.el.get_positions()
        for orb in orbitals:
            if os.path.isfile("basis_orbital_%i.hb" % (orb['index'])):
                print "found basis_orbital_%i" % orb['index']
                f = open("%s/basis_orbital_%i.hb" % (os.getcwd(), orb['index']),'r')
                values = pickle.load(f)
                f.close()
                assert nu.values(shape) == self.dim
            else:
                t1 = time()
                print "Calculating basis function %i/%i to grid..." % (orb['index']+1, self.calc.st.norb),
                atom = orb['atom']
                R = positions[atom]
                r_nl = orb['Rnl']
                l, m = self.transform(orb['orbital'])
                y_lm = Ylm(l, m)
                values = nu.zeros(self.dim, dtype=nu.complex)
                for i, x in enumerate(self.grid_points[0]):
                    for j, y in enumerate(self.grid_points[1]):
                        for k, z in enumerate(self.grid_points[2]):
                            vec = nu.array((x,y,z))-nu.array(R)
                            r, theta, phi = self.to_spherical_coordinates(vec)
                            values[i,j,k] = r_nl(r)*y_lm((theta, phi))
                print "in %i seconds." % (time() - t1)
                if self.store_grids:
                    f = open("%s/basis_orbital_%i.hb" % (os.getcwd(), orb['index']),'w')
                    pickle.dump(values, f)
                    f.close()
            self.basis_functions[orb['index']] = values


    def analyse_states(self):
        """ Performs the angular momentum analysis on all states. """
        for n in range(self.calc.st.norb):
            wf_coefficients = self.calc.st.wf[:,n]
            state_grid = nu.zeros(self.dim, dtype=nu.complex)
            for wf_coef, orb in zip(wf_coefficients, self.calc.el.orbitals()):
                state_grid += wf_coef * self.basis_functions[orb['index']]
            state_grid_squared = state_grid.conjugate() * state_grid
            self.norms[n] = nu.sum(state_grid_squared) * self.dV
            self.weights[n] = nu.sum(state_grid_squared * nu.where(self.shell_index_grid != 0, 1, 0)) * self.dV
            t1 = time()
            print "Analysing state %i..." % (n+1),
            for l in self.l_array:
                self.c_nl[n,l] = self.analyse_state(state_grid, l)
            print "in %i seconds." % (time() -  t1)
        for n in range(self.calc.st.norb):
            self.c_nl[n,:] = self.c_nl[n,:]/float(nu.sum(self.c_nl[n,:]))
        f = open("ylm_expansion_coefficients.hb",'w')
        pickle.dump(self.c_nl, f)
        f.close()


    def analyse_state(self, state_grid, l):
        """ Performs the angular momentum analysis with respect to
            angular momentum l to the state n. """
        c = 0.0
        for m in range(-l,l+1):
            ylm = self.ylms["%i,%i" % (l, m)]
            # The integration
            solid_angle_integrals = nu.zeros((len(self.shells.keys())), dtype=nu.float64)
            for i in self.shells.keys():
                # the shell grid
                shell_grid = nu.where(self.shell_index_grid == i, 1, 0)
                # the number of boxes in the i:th shell
                N_shell = self.shells[i]
                # the integration over the solid angle
                phi_nlm = nu.sum(shell_grid * ylm.conjugate() * state_grid)
                c += phi_nlm.conjugate() * phi_nlm * 4*pi*self.dV / N_shell
        return c


    def greetings(self):
        print "\n*** Starting the angular momentum analysis. ***"
        print "The grid contains %i x %i x %i grid points" % tuple(self.dim)
        print "The center of the expansion (in Ang): %0.2f, %0.2f, %0.2f" % tuple(self.origin * Bohr)
        print "The radius of the expansion (in Ang): %0.2f" % (self.R_0 * Bohr)
        print "The analysis is performed on angular momentums:",
        for l in self.l_array:
            print self.letters[l],
        print "\n"


    def write_to_file(self):
        e = self.calc.st.get_eigenvalues() * Hartree
        occ = self.calc.st.get_occupations()
        w = self.weights
        norms = self.norms
        f = open(self.filename, 'w')
        print >> f, "# The center of the expansion (in Ang): %0.2f, %0.2f, %0.2f" % tuple(self.origin * Bohr)
        print >> f, "# The radius of the expansion (in Ang): %0.2f" % (self.R_0 * Bohr)
        print >> f, "# The shell thickness (in Ang): %0.2f" % (self.a * Bohr)
        print >> f, "#state  energy(eV)   norm   weight     occ",
        for l in self.l_array:
            print >> f, "%6s" % self.letters[l],
        print >> f, ""
        for n in range(self.calc.st.norb):
            print >> f, "%5i %12.4f %7.4f %7.4f %7.4f" % (n, e[n], norms[n], w[n]/norms[n], occ[n]),
            for l in self.l_array:
                print >> f, "%6.3f" % self.c_nl[n,l],
            print >> f, ""
        f.close()


    def gaussian_peak(self, x, x0, width):
        return exp( - (x-x0)**2 / (4*width**2) )


    def make_plot(self, width=0.1, filename=None, bands=None):
        import pylab
        colors = ['#FFFF00','#FF0000','#5FD300','#2758D3',
                  '#058C00','#E1AB18','#50E1D0']

        if bands == None:
            bands = self.calc.st.norb
        if filename != None:
            f = open(filename,'r')
            self.c_nl = pickle.load(f)
            f.close()
            print "Loaded expansion coefficients from file."
        c_nl = self.c_nl
        e = self.calc.st.get_eigenvalues()[:bands]*Hartree
        e_min, e_max = min(e), max(e)
        empty = 0.1*(e_max - e_min)
        x_grid = nu.linspace(e_min-empty, e_max+empty, 400)
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
        pylab.show()


    def run(self):
        self.greetings()
        self.mark_grids()
        self.ylms_to_grid()
        self.basis_functions_to_grid()
        self.analyse_states()
        self.write_to_file()
        self.make_plot()


if __name__ == '__main__':
     import sys
     from ase import *
     from hotbit import Calculator
     if sys.argv[1] == 'H2':
         h2 = Atoms('H2', ((0,0,0),(1,0,0)))
         h2.center(vacuum=5)
         h2.set_calculator(Calculator(SCC=True))
         h2.get_potential_energy()
         JA = JelliumAnalysis(h2, maxl=1, R_0=3, a=0.4)
         JA.run()

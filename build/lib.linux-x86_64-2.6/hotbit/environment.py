import numpy as nu
from ase.units import Hartree
from weakref import proxy
hbar=0.02342178


class LinearlyPolarizedLaser:
    """ Class for laser potential. """

    def __init__(self,energy,flux,polarization,phase=0.0):
        """ Laser parameters.

        electrostatic potential=-E0*sin(omega*t+phase)*dot(polarization,r)

        Parameters:
        -----------
        omega: laser energy (in eV)
        flux: electric field flux  [ E_0=5.3E-9*sqrt(flux) ]
              flux ~ 10E8...10E16, ~10E12='normal' laser?
        polarization: direction for static polarization (3-array)
        phase: phase for laser pulse (do not start from zero)
        """
        self.omega=(energy/Hartree)/hbar # hbar*omega=energy
        self.E0=nu.sqrt(flux)*5.33802445585E-09
        self.pol=polarization/nu.linalg.norm(polarization)
        self.phase=phase



    def __call__(self,r,t):
        """ Return the electrostatic potential.

        Parameters:
        -----------
        r: position in Bohrs
        t: time in atomic units (~fs)
        """
        return -self.E0*nu.sin(self.omega*t+self.phase)*nu.dot(r,self.pol)




class Environment:
    def __init__(self,calc):
        self.t=0.0
        self.phis=[]
        self.calc=proxy(calc)

    def __del__(self):
        pass

    def propagate_time(self,dt):
        """ Increase time by dt (dt in atomic units) """
        self.t+=dt


    def add_phi(self,phi):
        """ Add external electrostatic potential function in atomic units.

        phi=phi(r,t) is any function, where r is position (Bohrs)
        and t time (atomic units, ~fs)
        """
        self.phis.append(phi)


    def phi(self,r):
        """ Return external electrostatic potential.

        Current internal environment time is used.

        Parameters:
        -----------
        r: position (Bohrs) or atom index;
           index if r is an integer
        """
        if isinstance(r,int):
            r=self.calc.el.nvector(r)
        pot=0.0
        for f in self.phis:
            pot+=f(r,self.t)
        return pot






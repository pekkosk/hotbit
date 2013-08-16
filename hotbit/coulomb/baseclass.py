"""
Coulomb base class.
"""


class Coulomb:
    def __init__(self):
        pass

    def get_potential(self, a=None):
        """
        Return the electrostatic potential for each atom.
        """
        raise NotImplementedError()

    def get_potential_and_field(self, a=None):
        """
        Return the both, the electrostatic potential and the field for each
        atom.
        """
        raise NotImplementedError()

    def get_gamma(self, a=None):
        """
        Return the gamma correlation matrix, i.e. phi(i) = gamma(i, j)*q(j).
        """
        raise NotImplementedError()

    def get_potential_energy(self, a=None):
        """
        Return the Coulomb energy.
        """
        raise NotImplementedError()

    def get_forces(self, a=None):
        """
        Return forces.
        """
        raise NotImplementedError()
    
    

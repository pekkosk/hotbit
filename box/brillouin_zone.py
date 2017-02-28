"""
Contains data for high symmetry points in different structures.
"""

import numpy as np

import ase.dft.kpoints as kpoints
from functools import reduce

_Gamma   = [ 0.00, 0.00, 0.00 ]
_X       = [ 1.00, 0.00, 0.00 ]
_W       = [ 1.00, 0.50, 0.00 ]
_K       = [ 0.75, 0.75, 0.00 ]
_L       = [ 0.50, 0.50, 0.50 ]

_fcc_a1  = [   0., 1./2, 1./2 ]
_fcc_a2  = [ 1./2,   0., 1./2 ]
_fcc_a3  = [ 1./2, 1./2,   0. ]
_fcc_A   = np.array( [ _fcc_a1, _fcc_a2, _fcc_a3 ] )

symmetry_points = {
    'fcc': {
        'Gamma': np.dot(_fcc_A, _Gamma),
        'X':     np.dot(_fcc_A, _X),
        'W':     np.dot(_fcc_A, _W),
        'K':     np.dot(_fcc_A, _K),
        'L':     np.dot(_fcc_A, _L)
        }
    }

band_structures = {
    'fcc': [ 'Gamma', 'X', 'W', 'K', 'Gamma', 'L', 'W' ]
    }


def get_bandpath(crystal, cell, npoints=50):
    """
    Return a list of k-points sampling the band-structure for the respective
    crystal.
    """
    bs = band_structures[crystal]
    symp = symmetry_points[crystal]

    points = reduce(
        lambda y1, y2: y1+y2,
        [ [ symp[x1], symp[x2] ] for x1, x2 in zip(bs[:-1], bs[1:]) ]
        )

    return kpoints.get_bandpath(points, cell, npoints)

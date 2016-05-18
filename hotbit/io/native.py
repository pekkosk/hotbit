"""
Native Hotbit file format.
"""
from copy import copy

import numpy as np

from box import mix
from box.data import data as box_data
from box.interpolation import Function



def read_element_from_elm(fileobj, symbol):
    """
    Read element data from Hotbit-style .elm file.

    Parameters:
    -----------
    fileobj:   filename of file-object to read from
    symbol:    chemical symbol of the element
    """

    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    symb = mix.find_value(fileobj, 'symbol').strip()
    assert symb == symbol

    data = copy(box_data[symbol])

    data['U'] = float(mix.find_value(fileobj, 'U'))
    data['FWHM'] = float(mix.find_value(fileobj, 'FWHM'))
    data['epsilon'] = { }
    for orbital in data['valence_orbitals']:
        data['epsilon'][orbital] = float(
            mix.find_value(fileobj, 'epsilon_%s' % orbital)
            )
    data['comment'] = mix.find_value(fileobj, 'comment',
                                     fmt='strings',
                                     default=['no comment'])

    functions = _read_functions(fileobj, data['valence_orbitals'])

    energies=[]
    for orbital in data['valence_orbitals']:
        eps = data['epsilon'][orbital]
        if   's' in orbital: n=1
        elif 'p' in orbital: n=3
        elif 'd' in orbital: n=5
        energies.extend( [eps]*n )
    data['onsite_energies'] = energies
    data['nr_basis_orbitals'] = len(energies)
    data['valence_energies'] = np.array(energies, dtype=float)

    # vdW correction
    data['C6'] = mix.find_value(fileobj, 'C6', default=-1.0)
    data['p']  = mix.find_value(fileobj, 'p',  default=-1.0)
    data['R0'] = mix.find_value(fileobj, 'R0', default=-1.0)

    return data, functions



def _read_functions(fileobj, valence_orbitals):
    """
    Read radial wave functions (R=u/r), self-consistent potentials,
    confinements, etc. from given file.
    """

    default = np.array([[1,0],[2,0],[3,0],[4,0]])
    functions = {
        'unl': {},
        'Rnl': {}
        }

    for nl in valence_orbitals:
        m = mix.find_value(fileobj, 'u_%s' % nl, fmt='matrix', default=default)
        functions['unl'][nl]=Function('spline', m[:,0], m[:,1])
        functions['Rnl'][nl]=Function('spline', m[:,0], m[:,1]/m[:,0])

    m = mix.find_value(fileobj, 'effective_potential',
                        fmt='matrix', default=default)
    functions['effective_potential'] = Function('spline', m[:,0], m[:,1])
    m = mix.find_value(fileobj, 'confinement_potential',
                       fmt='matrix', default=default)
    functions['confinement_potential'] = Function('spline', m[:,0], m[:,1])

    return functions



def read_HS_from_par(fileobj, symboli, symbolj):
    """
    Read Hamiltonian and overlap data from Hotbit-style .par file.

    Parameters:
    -----------
    fileobj:   filename of file-object to read from
    symboli:   chemical symbol of the first element
    symbolj:   chemical symbol of the second element
    """

    if mix.find_value(fileobj, 'X_X_table', fmt='test'):
        table = mix.find_value(fileobj,
                               'X_X_table' % ( symboli, symbolj ),
                               fmt='matrix')
    else:
        table = mix.find_value(fileobj,
                               '%s_%s_table' % ( symboli, symbolj ),
                               fmt='matrix')

    return table[:, 0], table[:, 1:]



def read_repulsion_from_par(fileobj):
    """
    Read repulsion data from Hotbit-style .par file.

    Parameters:
    -----------
    fileobj:   filename of file-object to read from
    """

    try:
        v = mix.find_value(fileobj, 'repulsion', fmt='matrix')
    except:
        v = np.array([[0,0],[1,0],[2,0],[3,0]])

    return v[:, 0], v[:, 1]

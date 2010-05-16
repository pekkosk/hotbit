"""
DFTB file format.
See: http://www.dftb.org/
"""

import numpy as np

from hotbit.io.fortran import fortran_readline

def read_HS_skf(fileobj, si, sj):
    """
    Read Hamitonian and overlap data from DFTB-style .skf file.

    Parameters:
    -----------
    fileobj:   filename of file-object to read from
    si:        chemical symbol of the first element
    sj:        chemical symbol of the second element
    """

    if isinstance(fileobj, str):
        fileobj = open(fileobj)
    
    # Homoatomic interactions also contain element data
    if si == sj:
        # First line contains grid spacing and number of grid points
        dx, n, dummy = fortran_readline(fileobj)
        n = int(n)

        # Contains self-energies, spin-polarization energies, Hubbard-U, ...
        l = fileobj.readline()
        # Don't know what this is for
        l = fileobj.readline()
    else:
        # First line contains grid spacing and number of grid points
        dx, n = fortran_readline(fileobj)
        n = int(n)

    x = dx*np.arange(0, n)

    HS = [ ]
    for i in range(n):
        HS += [ [ x[i] ] + fortran_readline(fileobj) ]
    HS = np.array(HS)

    return HS


def read_rep_skf(fileobj, rep_x0=0.1, rep_dx=0.005):
    """
    Read repulsion from DFTB-style .skf file.
    The repulsion in the .skf-file consists of an exponential and a spline
    part. Both will be converted to a single table.

    Parameters:
    -----------
    fileobj:   filename of file-object to read from
    rep_x0:    minimum value of the repulsion table
    rep_dx:    step size for discretization
    """

    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    l = fileobj.readline()
    while l and l.strip() != 'Spline':
        l = fileobj.readline()

    if not l:
        raise RuntimeError('Could not find "Spline" keyword when reading '
                           'repulsion.')

    n, cutoff = fortran_readline(fileobj)
    n = int(n)

    c1, c2, c3 = fortran_readline(fileobj)

    x1, x2, splc1, splc2, splc3, splc4 = fortran_readline(fileobj)

    x = np.linspace(rep_x0, cutoff, int((cutoff-rep_x0)/rep_dx)+1)
    i0 = np.searchsorted(x, x1)

    y = np.zeros(len(x))
    y[:i0] = c3 + np.exp(c2-c1*x[:i0])

    for j in range(n-1):
        if j > 0:
            last_x2 = x2
            x1, x2, splc1, splc2, splc3, splc4 = fortran_readline(fileobj)
            assert x1 == last_x2
        i1 = np.searchsorted(x, x2)
        y[i0:i1] = \
            splc1 + \
            splc2 * (x[i0:i1]-x1) + \
            splc3 * (x[i0:i1]-x1)**2 + \
            splc4 * (x[i0:i1]-x1)**3
        i0 = i1

    # The last entry is a fifth-order polynomial
    last_x2 = x2
    x1, x2, splc1, splc2, splc3, splc4, splc5, splc6 = fortran_readline(fileobj)
    assert x1 == last_x2

    i1 = np.searchsorted(x, x2)
    y[i0:i1] = \
        splc1 + \
        splc2 * (x[i0:i1]-x1) + \
        splc3 * (x[i0:i1]-x1)**2 + \
        splc4 * (x[i0:i1]-x1)**3 + \
        splc5 * (x[i0:i1]-x1)**4 + \
        splc6 * (x[i0:i1]-x1)**5

    return x, y

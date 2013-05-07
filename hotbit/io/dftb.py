"""
DFTB file format.
See: http://www.dftb.org/
"""
from math import log, pi, sqrt
from copy import copy

import numpy as np

from box.data import data as box_data

from hotbit.io.fortran import fortran_readline



def read_element_from_skf(fileobj, symbol):
    """
    Read element data from DFTB-style .skf file.

    Parameters:
    -----------
    fileobj:   filename of file-object to read from
    symbol:    chemical symbol of the element
    """

    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    # First line contains grid spacing and number of grid points, ignore
    fileobj.readline()

    # Self-energies, spin-polarization energies, Hubbard-U
    eself = [ 0.0 ]*3
    U = [ 0.0 ]*3
    q = [ 0.0 ]*3
    eself[0], eself[1], eself[2], espin, \
        U[0], U[1], U[2], q[0], q[1], q[2] = fortran_readline(fileobj)

    # Initialize data from database
    data = copy(box_data[symbol])

    data['epsilon'] = { }
    for orbital in data['valence_orbitals']:
        if 's' in orbital:
            data['epsilon'][orbital] = eself[2]
        elif 'p' in orbital:
            data['epsilon'][orbital] = eself[1]
        elif 'd' in orbital:
            data['epsilon'][orbital] = eself[0]

    # Apparently, only the last U is used
    data['U'] = U[2]
    data['FWHM'] = sqrt(8*log(2.0)/pi)/data['U']

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

    data['comment'] = None

    # .skf files do not contain basis information
    return data, None


def read_HS_from_skf(fileobj, symboli, symbolj):
    """
    Read Hamitonian and overlap data from DFTB-style .skf file.

    Parameters:
    -----------
    fileobj:   filename of file-object to read from
    symboli:   chemical symbol of the first element
    symbolj:   chemical symbol of the second element
    """

    if isinstance(fileobj, str):
        fileobj = open(fileobj)
    
    # Homoatomic interactions also contain element data
    if symboli == symbolj:
        # First line contains grid spacing and number of grid points
        l = fortran_readline(fileobj)
        dx = l[0]
        n = l[1]
        n = int(n)

        # Contains self-energies, spin-polarization energies, Hubbard-U, ...
        l = fileobj.readline()
    else:
        # First line contains grid spacing and number of grid points
        dx, n = fortran_readline(fileobj)
        n = int(n)

    # The n information is sometimes wrong, better count while reading
    #x = dx*np.arange(0, n)

    HS = [ ]
    l = fileobj.readline()
    while l and l.strip() != 'Spline':
        if l.strip() == 'Spline':
            if i != n-1:
                raise RuntimeError('Spline keyword encountered reading tables '
                                   'for %s-%s before reaching the end of the '
                                   'HS table.' % ( symboli, symbolj ))
        else:
            HS += [ fortran_readline(l) ]

        l = fileobj.readline()

    # don't care if not spline...
    #if not l:
    #    raise RuntimeError('Premature end-of-file: Keyword "Spline" not found '
    #                       'for %s-%s.' % ( symboli, symbolj ))
    
    HS = np.array(HS)
    x = dx*np.arange(0, HS.shape[0])

    #return x[0:HS.shape[0]], np.array(HS)
    return x, np.array(HS)



def read_repulsion_from_skf(fileobj, rep_x0=0.1, rep_dx=0.005):
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
        # no spline
        fileobj.seek(0)
        dx, n = fortran_readline(fileobj) #ignore
        items = fortran_readline(fileobj)
        if len(items)<20: # in homonuclear tables one element info came before repulsion
            items = fortran_readline(fileobj) 
        cutoff = items[9]
        c = [0,0] + items[1:9]         
        x = np.linspace(rep_x0, cutoff, int((cutoff-rep_x0)/rep_dx)+1)
        y = np.zeros_like(x)
        for i in range(2,10):
            y = y + c[i]*(cutoff-x)**i
        
        #raise RuntimeError('Could not find "Spline" keyword when reading '
        #                   'repulsion.')

    else:
        n, cutoff = fortran_readline(fileobj)
        n = int(n)
    
        # Coefficients for the exponential
        c1, c2, c3 = fortran_readline(fileobj)
    
        x1, x2, splc1, splc2, splc3, splc4 = fortran_readline(fileobj)
    
        x = np.linspace(rep_x0, cutoff, int((cutoff-rep_x0)/rep_dx)+1)
        y = c3 + np.exp(c2-c1*x)
    
        i0 = np.searchsorted(x, x1)+1
        for j in range(n-1):
            if j > 0:
                last_x2 = x2
                x1, x2, splc1, splc2, splc3, splc4 = fortran_readline(fileobj)
                assert x1 == last_x2
            i1 = np.searchsorted(x, x2)+1
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
    
        i1 = np.searchsorted(x, x2)+1
        y[i0:i1] = \
            splc1 + \
            splc2 * (x[i0:i1]-x1) + \
            splc3 * (x[i0:i1]-x1)**2 + \
            splc4 * (x[i0:i1]-x1)**3 + \
            splc5 * (x[i0:i1]-x1)**4 + \
            splc6 * (x[i0:i1]-x1)**5

    return x, y

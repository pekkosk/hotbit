"""
Construct neighbor lists and return distances and distance vectors.
"""

# Copyright (C) 2010 NSC Jyvaskyla, Fh-IWM
# Please see the accompanying LICENSE file for further information.

import numpy as np

# diag_indices_from was introduced in numpy 1.4.0
if hasattr(np, 'diag_indices_from'):
    diag_indices_from = np.diag_indices_from
else:
    def diag_indices_from(m):
        i = [ ] 
        for n in m.shape:
            i += [ np.arange(n, dtype=int) ]
        return tuple(i)


def n_from_ranges(s, n, m=None):
    r = [ ]

    if type(n) == int or type(n) == float:
        n = [n]*3
    if m is not None:
        if type(m) == int or type(m) == float:
            m = [m]*3
    else:
        m = [ i for i in n ]
        n = [ -i for i in n ]

    for cn, cm, ( i1, i2 ) in zip(n, m, s):
        if i1 == -np.Inf:
            j1 = cn
        else:
            j1 = int(i1)
        if i2 == np.Inf:
            j2 = cm+1
        else:
            j2 = int(i2)+1
        r += [ ( j1, j2 ) ]

    return r


def sappend(x, y):
    if x is None:
        return y
    elif len(x) == 0:
        return y
    else:
        return np.append(x, y, axis=0)


def get_neighbors(a, cutoff=None):
    """
    Given a Hotbit atoms object, return a list of neighbors that are within a
    certain cutoff. Neighbors are returned for each atom, along with the
    distance and the normal vector pointing to the neighbor.

    Parameters:
    -----------
    a:        Hotbit atoms object
    cutoff:   Return neighbors up to this cutoff. If omitted, all neighbors are
              returned.

    Return tuple:
    -------------
    i:   Indices of the primary atom
    j:   Indices of the neighboring atom
    d:   Distances
    n:   Normal vectors
    """
    # FIXME!!! Brute force scan in the neighboring 5 cells in all directions.
    # This might not be enough, i.e. for chiral nanotubes. Maybe the container
    # itself should contain a function to return nearest neighbor periodic
    # cells.
    # For large systems, some kind of binning scheme should implemented.
    sym_ranges  = a.get_symmetry_operation_ranges()
    if cutoff is None and not a.is_cluster():
        raise RuntimeError("Please specify a cutoff when searching for "
                           "neighbors in a periodic system.")

    if cutoff is None:
        n1, n2, n3  = n_from_ranges(sym_ranges, np.Inf)
    else:
        n1, n2, n3  = n_from_ranges(sym_ranges, 5)

    nat = len(a)
    r = a.get_positions()

    # Some reference coordinate, the origin should be fine
    r0_v = np.zeros(3, dtype=float)

    jl = np.tile(np.arange(nat, dtype=int), (nat, 1))
    il = jl.transpose()

    i = None
    j = None
    d = None
    n = None

    # Contribution of neighboring boxes
    for x1 in range(*n1):
        for x2 in range(*n2):
            for x3 in range(*n3):
                # construct a matrix with distances
                r1      = a.transform(r0_v, [x1, x2, x3])
                T       = a.rotation([x1, x2, x3])
                    
                rT      = np.dot(r-r0_v, np.transpose(T))

                dr      = r.reshape(nat, 1, 3) - (r1+rT).reshape(1, nat, 3)
                abs_dr  = np.sqrt(np.sum(dr*dr, axis=2))
                if cutoff is None:
                    mask  = np.ones([nat, nat], dtype=bool)
                else:
                    mask  = abs_dr < cutoff

                # Do not return self-interactions
                if x1 == 0 and x2 == 0 and x3 == 0:
                    mask[diag_indices_from(mask)] = False

                if np.any(mask):
                    i       = sappend(i, il[mask])
                    j       = sappend(j, jl[mask])

                    abs_dr  = abs_dr[mask]
                    dr      = dr[mask, :]/abs_dr.reshape(-1, 1)

                    d       = sappend(d, abs_dr)
                    n       = sappend(n, dr)

    return i, j, d, n

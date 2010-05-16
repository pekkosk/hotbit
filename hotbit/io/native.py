"""
Native Hotbit file format.
"""

from box import mix


def read_HS_par(fileobj, si, sj):
    """
    Read Hamiltonian and overlap data from Hotbit-style .par file.
    """
    if mix.find_value(fileobj, 'X_X_table', fmt='test'):
        table = mix.find_value(fileobj, 'X_X_table' % (si, sj), fmt='matrix')
    else:
        table = mix.find_value(fileobj, '%s_%s_table' % (si, sj), fmt='matrix')
    
    return table



def read_rep_par(fileobj):
    """
    Read repulsion data from Hotbit-style .par file.
    """
    try:
        v = mix.find_value(fileobj, 'repulsion', fmt='matrix')
    except:
        v = nu.array([[0,0],[1,0],[2,0],[3,0]])
    
    return v[:, 0], v[:, 1]

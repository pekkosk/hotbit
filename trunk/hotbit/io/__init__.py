"""
Input module. Contains functions to read element, Slater-Koster and repulsion
data.
"""

def read_element(filename, symbol, format=None):
    """
    Read element data from files.

    Parameters:
    -----------
    fileobj:   filename of file-object to read from
    symbol:    chemical symbol of the element
    """

    if format is None:
        format = filetype(filename)

    if format == 'elm':
        from hotbit.io.native import read_element_from_elm
        return read_element_from_elm(filename, symbol)

    if format == 'skf':
        from hotbit.io.dftb import read_element_from_skf
        return read_element_from_skf(filename, symbol)

    raise RuntimeError('File format "'+format+'" not recognized!')



def read_HS(filename, symboli, symbolj, format=None):
    """
    Read Slater-Koster tables from files.

    Parameters:
    -----------
    fileobj:   filename of file-object to read from
    symboli:   chemical symbol of the first element
    symbolj:   chemical symbol of the second element
    """

    if format is None:
        format = filetype(filename)

    if format == 'par':
        from hotbit.io.native import read_HS_from_par
        return read_HS_from_par(filename, symboli, symbolj)

    if format == 'skf':
        from hotbit.io.dftb import read_HS_from_skf
        return read_HS_from_skf(filename, symboli, symbolj)

    raise RuntimeError('File format "'+format+'" not recognized!')



def read_repulsion(filename, format=None):
    """
    Read Slater-Koster tables from files.

    Parameters:
    -----------
    fileobj:   filename of file-object to read from
    """

    if format is None:
        format = filetype(filename)

    if format == 'par':
        from hotbit.io.native import read_repulsion_from_par
        return read_repulsion_from_par(filename)

    if format == 'skf':
        from hotbit.io.dftb import read_repulsion_from_skf
        return read_repulsion_from_skf(filename)

    raise RuntimeError('File format "'+format+'" not recognized!')



def filetype(filename):
    if filename.lower().endswith('.elm'):
        return 'elm'

    if filename.lower().endswith('.par'):
        return 'par'

    if filename.lower().endswith('.skf') or filename.lower.endswith('.spl'):
        return 'skf'

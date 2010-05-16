

def read_HS(filename, si, sj, format=None):
    """
    Read Slater-Koster tables from files.
    """

    if format is None:
        format = filetype(filename)

    if format == 'par':
        from hotbit.io.native import read_HS_par
        return read_HS_par(filename, si, sj)

    if format == 'skf':
        from hotbit.io.dftb import read_HS_skf
        return read_HS_skf(filename, si, sj)

    raise RuntimeError('File format "'+format+'" not recognized!')


def read_rep(filename, format=None):
    """
    Read Slater-Koster tables from files.
    """

    if format is None:
        format = filetype(filename)

    if format == 'par':
        from hotbit.io.native import read_rep_par
        return read_rep_par(filename)

    if format == 'skf':
        from hotbit.io.dftb import read_rep_skf
        return read_rep_skf(filename)

    raise RuntimeError('File format "'+format+'" not recognized!')


def filetype(filename):
    if filename.lower().endswith('.par'):
        return 'par'

    if filename.lower().endswith('.skf') or filename.lower.endswith('.spl'):
        return 'skf'

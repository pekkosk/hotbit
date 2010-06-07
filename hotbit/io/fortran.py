"""
Input of Fortran formatted text files.
"""

def fortran_readline(f, dtype=float):
    """
    Read a line from a fortran text file and convert to a list.

    In particular, statements like '5*20.0' are resolved to a list containing
    5 times the value 20.0.
    """

    if isinstance(f, str):
        l = f
    else:
        l = f.readline()
    s = l.replace('d', 'e').replace(',', ' ').split()

    r = [ ]
    for e in s:
        i = e.find('*')

        if i == -1:
            r += [ dtype(e) ]
        else:
            n = int(e[:i])
            v = dtype(e[i+1:])

            r += [ v ]*n
    
    return r

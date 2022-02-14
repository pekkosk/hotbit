#!/usr/bin/python

"""
This script transforms Hotbit Parameters for specified Atoms in
Hotbit-Parameter files (*.par, *.elm) to DFTB (*.skf) files.
argv[1], argv[2] = atomic symbol (i.e. H H, C H, ..)

    usage: hb2skl C C

RJ Maurer & GS Michelitsch, Technische Universitaet Muenchen, 14/03/2014
"""

from sys import argv
from hotbit.io import native

def generate_DFTBplus_repulsion(filename, DEBUG=False):
    """
    Transforms Hotbit repulsion to DFTB+ format.
    This function generates output to be added to an .skf
    file. It outputs a string that has to be inserted
    in the right position.
    RJ Maurer, 14/03/2014
    """

    import numpy as np
    from scipy.interpolate import splrep, splev
    from scipy.optimize import curve_fit

    tmp = filename
    line = 0
    #read until the repulsion begins
    while not 'repulsion=' in tmp[line]:
        line += 1
    else:
        #skip commentary
        line += 1
        x, y = [], []
        while line < len(tmp):
            try:
                is_line = '.' in tmp[line].strip('\n').split()[0]
            except:
                break
            if not is_line:
                break
            data = tmp[line].strip('\n').split()
            x.append(float(data[0]))
            y.append(float(data[1]))
            line += 1

    ##Repulsion is defined by an initial exponential part
    #and a third order spline with weird knots with
    #a 5th order spline at the end

    #repulsion function is stored in x and y
    def expfit(x, a1, a2, a3):
        return np.exp(-a1 * x + a2) + a3

    def splinefit(t, c0, c1, c2, c3):
        #here t = (x-x0) directly
        func = c0 + c1*t + c2*t**2 + c3*t**3
        return func

    def splinefit2(t, c0, c1, c2, c3, c4, c5):
        # here t = (x-x0) directly
        func = c0 + c1 * t + c2 * t ** 2 + \
            c3 * t ** 3 + c4 * t ** 4 + \
            c5 * t ** 5
        return func

    #third order scipy spline
    SPLINE = splrep(x, y, s=0, k=3)
    if DEBUG:
        import matplotlib.pylab as pl
        pl.axhline(0,c='k',ls='--')
        pl.plot(x,y,lw=4,ls='--',c='r',alpha=0.2,label='original')

    #START WITH EXPONENTIAL FIT
    #find last x value below 0.50 Bohr
    s1 = x[0]
    x_cut = x[-1]
    for i in x:
        if i < 0.50:
            s2 = i

    print('Exponential function will be fitted in ' + \
        'the range {0} to {1}'.format(s1, s2))

    t = np.linspace(s1, s2, 10)
    data = splev(t, SPLINE)

    expParams, fitCovariances = curve_fit(expfit, t, data)
    print(' Exp. fit coefficients:\n', expParams)
    print(' Covariance matrix:\n', fitCovariances)
    repulsion_string = 'Spline\n'

    if DEBUG:
        import matplotlib.pylab as pl
        a1, a2, a3 = expParams
        X = np.linspace(s1, s2, 30)

        Y = [splev(x, SPLINE) for x in X]
        Y2 = [expfit(x, a1, a2, a3) for x in X]

        pl.plot(X, Y)
        pl.plot(X, Y2)

    #now we fit the cubic splines piecewise,except the last
    #hereby we split the range in pieces of length 0.10 Bohr,
    #except the last one with 0.4 Bohr
    deltax = 0.10
    final_dx = 0.30
    x_end = x_cut - final_dx
    N_dftbplus = int((x_end - s2) / deltax)
    x_end = s1 + deltax * N_dftbplus

    repulsion_string += '{0} {1:4.2f}\n'.format(N_dftbplus + 1, x_cut)
    a1, a2, a3 = expParams
    repulsion_string += '{0:18.15f} {1:18.15f} {2:18.15f}\n'.format(a1, a2, a3)

    start = s2
    end = s2 + deltax
    for knot in range(N_dftbplus):
        print('Fitting a spline basis function in the range between ' + \
            '{0} and {1}'.format(start, end))
        t = np.linspace(start, end, 6)
        data = splev(t, SPLINE)
        splParams, fitCovariances = curve_fit(splinefit, t-start, data)
        print(' Spline fit coefficients:\n', splParams)
        print(' Covariance matrix:\n', fitCovariances)

        c0, c1, c2, c3 = splParams
        repulsion_string += '{0:5.3f} {1:5.3f} {2:19.15f} {3:19.15f} \
{4:19.15f} {5:19.15f}\n'.format(start, end, c0, c1, c2, c3)

        if DEBUG:
            X = np.linspace(start, end, 20)
            Y = [splev(x, SPLINE) for x in X]
            Y2 = [splinefit(x-start, c0, c1, c2, c3) for x in X]
            pl.plot(X, Y)
            pl.plot(X, Y2)

        start = end
        end += deltax

    end = x_cut
    print('Now Fitting the final spline basis function in the range ' + \
          'between {0} and {1}'.format(start, end))

    t = np.linspace(start, end, 12)
    data = splev(t, SPLINE)
    splParams, fitCovariances = curve_fit(splinefit2,t-start,data)
    print(' Spline fit coefficients:\n', splParams)
    print(' Covariance matrix:\n', fitCovariances)

    c0, c1, c2, c3, c4, c5 = splParams
    repulsion_string += '{0:5.3f} {1:5.3f} {2:19.15f} {3:19.15f} {4:19.15f} \
{5:19.15f} {6:19.15f} {7:19.15f}\n'.format(start, end, c0, c1, c2, c3, c4, c5)

    if DEBUG:
        X = np.linspace(start, end, 20)
        Y = [splev(x, SPLINE) for x in X]
        Y2 = [splinefit2(x-start, c0, c1, c2, c3, c4, c5) for x in X]
        pl.plot(X, Y)
        pl.plot(X, Y2)
        pl.legend()
        pl.show()

    return repulsion_string




def generate_DFTBplus_header(elmdat, pardat, hetero):
    """
    Transforms Hotbit Metadata and Hamiltonian to DFTB+ format.
    This function generates output to be added to an .skf
    file. It outputs a string that has to be inserted
    in the right position.
    GS Michelitsch, 14/03/2014
    """
    # First Line
    # [gridDist]   [nGridPoints]   [?]
    # Uncommented non-documented 3rd parameter, it breaks Hotbit
    header = str("%.12f" % (pardat[0][1] - pardat[0][0])) + ', ' + \
        str(len(pardat[1]) + int(pardat[0][0] / (pardat[0][1] - pardat[0][0])))
    #+ ', ' + '0\n'

    header += '\n'

    if hetero:
        # Second Line
        # [Ed]   [Ep]   [Es]
        # (if present, otherwise 0.0)
        for i in range(3 - len(list(elmdat[0]['epsilon'].keys()))):
            header += '0.0 '
        for key in sorted(elmdat[0]['epsilon'].keys()):
            header += str("%.8f" % elmdat[0]['epsilon'][key]) + ' '

        # Second Line
        # [SPE]   [Ud]   [Up]   [Us]
        # (we don't have the SPE value in Hotbit: Spin Polarization Error of
        # early implementations in DFTB. Apparently only Us is evaluated,
        # therefore we have the same value for Ud Up Us)
        header += ', 0.0, ' + (str("%.6f" % elmdat[0]['U']) + ' ') * 3

        # Second Line
        # [fd]   [fp]   [fs]
        # Occupation numbers of valence orbitals
        for i in range(3 - len(elmdat[0]['valence_orbitals'])):
            header += '0.0 '
        for key in sorted(elmdat[0]['configuration'].keys()):
            if key in elmdat[0]['valence_orbitals']:
                header += str("%.1f" % elmdat[0]['configuration'][key]) + ' '

        header += '\n'

    # Third line
    # [mass]   [c2] - [c9]   [rcut]   [d1] - [d10]
    # Except for [mass], all parameters are set to 1.0
    header += str("%.2f" % elmdat[0]['mass']) + ', 19*1.0,\n'

    # Fourth Line
    # [Hxxx]
    # Integral table containing the DFTB Hamiltonian
    if pardat[0][0] != 0:
        for i in range(int(pardat[0][0] / (pardat[0][1] - pardat[0][0]))):
            header += str(len(pardat[1][0])) + '*1.0,\n'
    ct, theader = 0, ''
    for i in range(len(pardat[1])):
        for j in range(len(pardat[1][i])):
            if pardat[1][i][j] == 0:
                ct += 1
                theader = str(ct) + '*0.0 '
            else:
                ct = 0
                header += theader
                theader = ''
                header += '{0: 1.12e}  '.format(pardat[1][i][j])
        if theader != '':
            ct = 0
            header += theader
            theader = ''
            header += '{0: 1.12e}  '.format(pardat[1][i][j])
        header += '\n'

    return header


elmf = None

# Open Files
try:
    parf = open(argv[1] + '_' + argv[2] + '.par')
    elmf = open(argv[1] + '.elm')
except:
    parf = open(argv[2] + '_' + argv[1] + '.par')


# Read line-by-line (needed for repulsive term)
tmp = parf.readlines()


for carousel in [(argv[1], argv[2]), (argv[2], argv[1])]:

    # Open Outfile
    skff = open(carousel[0] + '-' + carousel[1] + '.skf', 'w')

    # Read relevant data (A-B and B-A is different!)
    elmdat = native.read_element_from_elm(elmf, argv[1])
    pardat = native.read_HS_from_par(parf, carousel[0], carousel[1])

    # Generate DFTBplus format string
    print('Reading in Hotbit Metadata and repulsive potential')
    header = generate_DFTBplus_header(elmdat,
                                      pardat,
                                      argv[1] == argv[2]) + \
        generate_DFTBplus_repulsion(tmp,
                                    DEBUG=False)

    skff.write(header)
    skff.close()

# Close files
parf.close()
elmf.close()

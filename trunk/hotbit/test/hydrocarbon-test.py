from ase import *
from hotbit import Hotbit
import numpy as nu
import sys

"""
This is a test script that gives insight how good are
the parametrizations for H-H, C-C and C-H. Reference
values are taken from

    D. Porezag et al., Physical Review B, volume 51, number 19, 1995

Same molecules are optimized using the default parametrizations
and the quantities are compared to the reference values.
"""


porezag_data = {
    'H2':{'HH':0.765},
    'CH':{'CH':1.138},
    'CH2':{'CH':1.134, 'HCH':98.6},
    'CH3':{'CH':1.114, 'HCH':116.8},
    'CH4':{'CH':1.116},
    'C2H2':{'CC':1.206, 'CH':1.099},
    'C2H4':{'CC':1.321, 'CH':1.113, 'CCH':116.3},
    'C2H6':{'CC':1.503, 'CH':1.119, 'HCH':108.0},
    'C3H4':{'C1C2':1.318, 'C2C3':1.509, 'C1H':1.109, 'HC1C2':148.4},
    'C3H6':{'CC':1.503, 'CH':1.114},
    'C4H10':{'C1C2':1.511, 'C2C3':1.520},
    'C6H6':{'CC':1.389, 'CH':1.114}
    }

reference_data = {
    'H2':{'HH':0.765},
    'CH':{'CH':1.152},
    'CH2':{'CH':1.135, 'HCH':99.1},
    'CH3':{'CH':1.093, 'HCH':120},
    'CH4':{'CH':1.101},
    'C2H2':{'CC':1.212, 'CH':1.078},
    'C2H4':{'CC':1.331, 'CH':1.098, 'CCH':116.4},
    'C2H6':{'CC':1.513, 'CH':1.105, 'HCH':107.2},
    'C3H4':{'C1C2':1.305, 'C2C3':1.510, 'C1H':1.091, 'HC1C2':149.5},
    'C3H6':{'CC':1.504, 'CH':1.095},
    'C4H10':{'C1C2':1.517, 'C2C3':1.532},
    'C6H6':{'CC':1.396, 'CH':1.095}
    }

def length(atoms, a, b):
    return nu.linalg.norm(atoms.positions[a] - atoms.positions[b])


def angle(atoms, a, b, c, unit='deg'):
    """ Return the angle a-b-c. """
    vec1 = atoms.positions[b] - atoms.positions[a]
    vec2 = atoms.positions[b] - atoms.positions[c]
    angle = nu.arccos(nu.dot(vec1, vec2)/(nu.linalg.norm(vec1)*nu.linalg.norm(vec2)))*180/nu.pi
    if unit == 'deg':
        k = 1.
    elif unit == 'rad':
        k = 180/nu.pi
    return k*angle


def optimize(atoms, fmax=0.01):
    atoms.center(vacuum=6)
    atoms.rattle()
    if atoms.get_calculator() == None:
        atoms.set_calculator(Hotbit())
    dyn = QuasiNewton(atoms, maxstep=0.01)
    dyn.run(fmax=fmax)


def H2(data):
    name = 'H2'
    data[name] = {}
    atoms = Atoms('H2',((0,    0,    0),
                        (0.76, 0,    0)))
    optimize(atoms)

    data[name]['HH'] = length(atoms, 0, 1)


def CH(data):
    name = 'CH'
    data[name] = {}
    atoms = Atoms('CH',((0,   0,    0),
                        (1.1, 0,    0)))
    optimize(atoms)

    data[name]['CH'] = length(atoms, 0, 1)


def CH2(data):
    name = 'CH2'
    data[name] = {}
    atoms = Atoms('CH2',
      ((7.348,      5.390,      5.640),
       (6.288,      5.745,      5.850),
       (8.038,      6.169,      6.098)))
    optimize(atoms)

    data[name]['CH'] = length(atoms, 0, 1)
    data[name]['HCH'] = angle(atoms, 1,0,2)


def CH3(data):
    name = 'CH3'
    data[name] = {}
    a = 1.1
    atoms = Atoms('CH3',
      ((6.550,      6.952,      6.000),
       (7.650,      6.952,      5.999),
       (6.001,      7.904,      5.999),
       (6.001,      6.000,      5.999)))
    optimize(atoms)

    data[name]['CH'] = length(atoms, 0, 1)
    data[name]['HCH'] = angle(atoms, 1, 0, 2)


def CH4(data):
    name = 'CH4'
    data[name] = {}
    atoms = Atoms('CH4',
      ((6.459,      6.786,      6.250),
       (6.422,      6.706,      7.388),
       (7.547,      6.808,      5.905),
       (5.935,      7.745,      5.918),
       (5.931,      5.886,      5.787)))
    optimize(atoms)
    data[name]['CH'] = length(atoms, 0, 1)


def C2H2(data):
    name = 'C2H2'
    data[name] = {}
    atoms = Atoms('C2H2',
      ((7.226,      5.993,      6.042),
       (8.382,      5.973,      6.040),
       (6.157,      6.012,      6.045),
       (9.451,      5.954,      6.037)))
    optimize(atoms)

    data[name]['CC'] = length(atoms, 0, 1)
    data[name]['CH'] = length(atoms, 0, 2)


def C2H4(data):
    name = 'C2H4'
    data[name] = {}
    atoms = Atoms('C2H4',
      ((6.469,      6.865,      6.009),
       (7.736,      6.864,      5.992),
       (5.902,      7.807,      6.017),
       (5.903,      5.922,      6.017),
       (8.303,      7.807,      5.984),
       (8.302,      5.921,      5.984)))
    optimize(atoms)

    data[name]['CC'] = length(atoms, 0, 1)
    data[name]['CH'] = length(atoms, 0, 2)
    data[name]['CCH'] = angle(atoms, 1, 0, 2)


def C2H6(data):
    name = 'C2H6'
    data[name] = {}
    atoms = Atoms('C2H6',
      ((6.434,      6.752,      6.343),
       (6.433,      6.751,      7.856),
       (7.506,      6.751,      5.987),
       (5.898,      7.680,      5.985),
       (5.898,      5.823,      5.987),
       (7.505,      6.751,      8.214),
       (5.896,      7.679,      8.213),
       (5.897,      5.822,      8.212)))
    optimize(atoms)


    data[name]['CC'] = length(atoms, 0, 1)
    data[name]['CH'] = length(atoms, 0, 3)
    data[name]['HCH'] = angle(atoms, 5, 1, 6)


def C3H4(data):
    name = 'C3H4'
    a=1.2
    data[name] = {}
    atoms = Atoms('C3H4',
      ((6.722,      6.495,      6.659),
       (7.979,      6.359,      6.670),
       (7.498,      7.780,      6.624),
       (5.748,      6.021,      6.667),
       (8.828,      5.687,      6.693),
       (7.555,      8.386,      7.548),
       (7.565,      8.327,      5.665)))
    optimize(atoms)

    data[name]['C1C2']  = length(atoms, 0, 1)
    data[name]['C2C3']  = length(atoms, 0, 2)
    data[name]['C1H']   = length(atoms, 0, 3)
    data[name]['HC1C2'] = angle(atoms, 3,0,1)


def C3H6(data):
    name = 'C3H6'
    data[name] = {}
    atoms = Atoms('C3H6',
      ((6.283,      6.485,      6.760),
       (7.603,      6.293,      6.078),
       (7.037,      7.670,      6.239),
       (6.267,      6.352,      7.851),
       (8.457,      6.032,      6.721),
       (7.518,      8.317,      6.988),
       (5.388,      6.181,      6.197),
       (7.579,      5.865,      5.065),
       (6.639,      8.150,      5.333)))
    optimize(atoms)

    data[name]['CC'] = length(atoms, 0, 1)
    data[name]['CH'] = length(atoms, 0, 3)


def C4H10(data):
    name = 'C4H10'
    data[name] = {}
    atoms = Atoms('C4H10',
      ((7.350,      7.576,      7.578),
       (8.621,      6.992,      7.003),
       (9.880,      7.612,      7.587),
      (11.151,      7.021,      7.018),
       (7.340,      8.686,      7.372),
       (7.338,      7.376,      8.689),
       (8.630,      5.900,      7.232),
       (8.619,      7.183,      5.903),
       (9.874,      8.703,      7.348),
       (9.880,      7.432,      8.689),
      (11.155,      5.912,      7.227),
      (11.168,      7.219,      5.907),
       (6.478,      7.068,      7.072),
      (12.023,      7.526,      7.527)))
    optimize(atoms)

    data[name]['C1C2'] = length(atoms, 0, 1)
    data[name]['C2C3'] = length(atoms, 1, 2)


def C6H6(data):
    name = 'C6H6'
    data[name] = {}
    a=1.2
    b=2.1
    atoms = Atoms('C6H6',
      ((9.444,      7.815,      6.065),
       (8.770,      8.981,      6.062),
       (7.423,      8.981,      6.061),
       (6.750,      7.814,      6.063),
       (7.424,      6.648,      6.062),
       (8.771,      6.649,      6.063),
      (10.533,      7.816,      6.069),
       (9.314,      9.925,      6.059),
       (6.879,      9.923,      6.057),
       (5.661,      7.814,      6.066),
       (6.879,      5.705,      6.063),
       (9.316,      5.706,      6.063)))
    optimize(atoms)

    data[name]['CC'] = length(atoms, 0, 1)
    data[name]['CH'] = length(atoms, 0, 6)


def geometry_test(data):
    C3H4(data)
    H2(data)
    CH(data)
    CH2(data)
    CH3(data)
    CH4(data)
    C2H2(data)
    C2H4(data)
    C2H6(data)
    C3H6(data)
    C4H10(data)
    C6H6(data)


order = ['H2','CH','CH2','CH3','CH4','C2H2','C2H4','C2H6','C3H4','C3H6','C4H10','C6H6']


def print_results(data):
    rel_error_current = 0.0
    rel_error_porezag = 0.0
    print ""
    print "Current: value of the quantity calculated with default parameters."
    print "Porezag: values reported in Phys. Rev. B, volume 51, number 19, 1995"
    print "LSD:     DFT_LSD values obtained from the same paper."
    print "d:       the relative error of the previous value compared to"
    print "         the LSD value (d=(x - x_LSD)/x_LSD)."
    print ""
    print "The two-atom variables are bond lengths (in Angtroms) and"
    print "three-atom variables angles (in degrees)."
    print ""
    print "%9s %9s %9s %9s %9s %9s %9s" % ('Molecule', 'Variable', 'LSD', 'Current', 'd', 'Porezag', 'd')
    N = 0
    for n, k in enumerate(order):
        for i, kv in enumerate(data[k].items()):
            N += 1
            key, value = kv
            if i == 0:
                molecule = k
            else:
                molecule = ''
            d = value/reference_data[k][key]-1
            rel_error_current += d**2
            print "%9s %9s %9.3f %9.3f %9.3f" % (molecule, key, reference_data[k][key], value, d),
            d = porezag_data[k][key]/reference_data[k][key]-1
            rel_error_porezag += d**2
            print "%9.3f %9.3f" % (porezag_data[k][key], d),
            print ''
    print " The RMS of relative errors:                %0.3f               %0.3f" % (nu.sqrt(rel_error_current/N), nu.sqrt(rel_error_porezag/N))


if __name__ == '__main__':
    data = {}
    geometry_test(data)
    print_results(data)

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
    atoms = Atoms('H2',((0,0,0),(0.76,0,0)))
    optimize(atoms)

    data[name]['HH'] = length(atoms, 0, 1)


def CH(data):
    name = 'CH'
    data[name] = {}
    atoms = Atoms('CH',((0,0,0),(1.1,0,0)))
    optimize(atoms)

    data[name]['CH'] = length(atoms, 0, 1)


def CH2(data):
    name = 'CH2'
    data[name] = {}
    atoms = Atoms('HCH',((-1,0,0),(0,0,0),(1,0,0)))
    optimize(atoms)

    data[name]['CH'] = length(atoms, 0, 1)
    data[name]['HCH'] = angle(atoms, 0,1,2)


def CH3(data):
    name = 'CH3'
    data[name] = {}
    a = 1.1
    atoms = Atoms('CH3', ((0,0,0),(a,0,0),(a*cos(2*pi/3), a*sin(2*pi/3),0),(a*cos(4*pi/3), a*sin(4*pi/3),0)))
    optimize(atoms)

    data[name]['CH'] = length(atoms, 0, 1)
    data[name]['HCH'] = angle(atoms, 1, 0, 2)


def CH4(data):
    name = 'CH4'
    data[name] = {}
    atoms = Atoms('CH4',((0,0,0),(0,0,1),(sin(2*pi/3), 0, cos(2*pi/3)),
           (cos(2*pi/3)*sin(2*pi/3), sin(2*pi/3)*sin(2*pi/3), cos(2*pi/3)),
           (cos(4*pi/3)*sin(2*pi/3), sin(4*pi/3)*sin(2*pi/3), cos(2*pi/3))))
    optimize(atoms)
    data[name]['CH'] = length(atoms, 0, 1)


def C2H2(data):
    name = 'C2H2'
    data[name] = {}
    atoms = Atoms('HCCH', ((0,0,0),(1.2,0,0),(2.4,0,0),(3.6,0,0)))
    optimize(atoms)

    data[name]['CC'] = length(atoms, 1, 2)
    data[name]['CH'] = length(atoms, 0, 1)


def C2H4(data):
    name = 'C2H4'
    data[name] = {}
    atoms = Atoms('H2C2H2',((cos(2*pi/3),sin(2*pi/3),0),
                            (cos(4*pi/3),sin(4*pi/3),0), (0,0,0),
                            (1.2,0,0),(1.2-cos(2*pi/3),sin(2*pi/3),0),
                            (1.2-cos(4*pi/3),sin(4*pi/3), 0)))
    optimize(atoms)

    data[name]['CC'] = length(atoms, 2, 3)
    data[name]['CH'] = length(atoms, 1, 2)
    data[name]['CCH'] = angle(atoms, 0, 2, 3)


def C2H6(data):
    name = 'C2H6'
    data[name] = {}
    atoms = Atoms('H3C2H3',((sin(2*pi/3), 0, cos(2*pi/3)),
                            (cos(2*pi/3)*sin(2*pi/3), sin(2*pi/3)*sin(2*pi/3), cos(2*pi/3)),
                            (cos(4*pi/3)*sin(2*pi/3), sin(4*pi/3)*sin(2*pi/3), cos(2*pi/3)),
                            (0, 0, 0),
                            (0, 0, 1.2),
                            (cos(0)*sin(2*pi/3), 0, 1.2-cos(2*pi/3)),
                            (cos(2*pi/3)*sin(2*pi/3), sin(2*pi/3)*sin(2*pi/3), 1.2-cos(2*pi/3)),
                            (cos(4*pi/3)*sin(2*pi/3), sin(4*pi/3)*sin(2*pi/3), 1.2-cos(2*pi/3))))
    optimize(atoms)


    data[name]['CC'] = length(atoms, 3, 4)
    data[name]['CH'] = length(atoms, 0, 3)
    data[name]['HCH'] = angle(atoms, 0, 3, 2)


def C3H4(data):
    name = 'C3H4'
    data[name] = {}
    atoms = Atoms('C3H4', ((0,0,0),(1.2,0,0),(0.6,0.8,0),
                           (-0.5,-0.5,0),(1.7,-0.5,0),
                           (0.6,1.3,0.5),(0.6,1.3,-0.5)))
    optimize(atoms)

    data[name]['C1C2']  = length(atoms, 0, 1)
    data[name]['C2C3']  = length(atoms, 1, 2)
    data[name]['C1H']   = length(atoms, 0, 3)
    data[name]['HC1C2'] = angle(atoms, 1,0,3)


def C3H6(data):
    name = 'C3H6'
    data[name] = {}
    atoms = Atoms('C3H6',((0,0,0),(1.2,0,0),(0.6,1,0),
                          (-0.5,-0.5,0.5),(1.7, -0.5, 0.5),(0.6, 1.5, 0.5),
                          (-0.5,-0.5,-0.5),(1.7,-0.5,-0.5),(0.6, 1.5,-0.5)))
    optimize(atoms)

    data[name]['CC'] = length(atoms, 0, 1)
    data[name]['CH'] = length(atoms, 0, 3)


def C4H10(data):
    name = 'C4H10'
    data[name] = {}
    atoms = Atoms('CH2CH2CH2CH2H2',((0,0,0),(0,1,0),(0,0,1),
                                    (1,0,0),(1,-1,0),(1,0,-1),
                                    (2,0,0),(2,1,0),(2,0,1),
                                    (3,0,0),(3,-1,0),(3,0,-1),
                                    (-1,0,0),(4,0,0)))
    atoms.positions = atoms.positions * 1.3
    optimize(atoms)

    data[name]['C1C2'] = length(atoms, 0, 3)
    data[name]['C2C3'] = length(atoms, 3, 6)


def C6H6(data):
    name = 'C6H6'
    data[name] = {}
    a=1.2
    b=2.1
    atoms = Atoms()
    for d in nu.arange(6) / 6. * 2*pi:
        atoms += Atom('C',(a*cos(d), a*sin(d), 0))
        atoms += Atom('H',(b*cos(d), b*sin(d), 0))
    optimize(atoms)

    data[name]['CC'] = length(atoms, 0, 2)
    data[name]['CH'] = length(atoms, 0, 1)


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

# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

from math import sqrt

import numpy as np
from ase.units import Hartree, Bohr
from weakref import proxy


def C6_from_p(p, nvel):
    N = 1.17 + 0.33*nvel
    return 0.75*sqrt(N*p)


def p_from_C6(C6, nvel):
    N = 1.17 + 0.33*nvel
    return (C6/0.75)**2/N


class vdWPairCorrection:

    def __init__(self, C6a, pa, nvela, R0a, C6b=None, pb=None, nvelb=None, R0b=None, d=3.0, N=7, M=4, Rmax=10.0):
        self.C6a   = C6a
        self.pa    = pa
        self.nvela = nvela
        self.R0a   = R0a
        if C6b is None and pb is None and nvelb is None:
            self.C6b   = C6a
            self.pb    = pa
            self.nvelb = nvela
        else:
            self.C6b   = C6b
            self.pb    = pb
            self.nvelb = nvelb
        if R0b is None:
            self.R0b = R0a
        else:
            self.R0b = R0b
        self.d = d
        self.N = N
        self.M = M
        self.Rmax = Rmax

        # Compute missing coefficients
        if self.C6a is None:
            self.C6a = C6_from_p(self.pa, self.nvela)
        elif self.pa is None:
            self.pa = p_from_C6(self.C6a, self.nvela)
        if self.C6b is None:
            self.C6b = C6_from_p(self.pb, self.nvelb)
        elif self.pb is None:
            self.pb = p_from_C6(self.C6b, self.nvelb)

        # Compute cross terms
        self.C6ab = 2 * self.C6a * self.C6b * self.pa * self.pb / (self.pa**2 * self.C6a + self.pb**2 * self.C6b)
        self.R0ab = (self.R0a**3 + self.R0b**3)/(self.R0a**2 + self.R0b**2)


    def __call__(self, r, der=0):
        if r is None:
            return self.Rmax

        h2 = self.d/(self.R0ab**self.N)
        h1 = np.exp(-h2*r**self.N)
        f  = (1.0-h1)**self.M
        h3 = self.C6ab/(r**6)

        if der == 0:
            return -f * h3
        elif der == 1:
            df = self.M*(1.0-h1)**(self.M-1) * h2*self.N*r**(self.N-1)*h1
            return ( 6 * f / r - df ) * h3


def setup_vdw(calc):
    #if calc.get('vdw'):
    #    raise NotImplementedError('van der Waals interactions are not yet implemented.')
    
    elms = len(calc.el.present)
    for i,s1 in enumerate(calc.el.present):
        for s2 in calc.el.present[i:]:      
            #
            # here vdw is simply the interaction
            # between elements s1 and s2
            #
            e1 = calc.el.elements[s1]
            e2 = calc.el.elements[s2]
            vdW = vdWPairCorrection(e1.get_C6(), e1.get_p(), e1.get_valence_number(), e1.get_R0(),
                                    e2.get_C6(), e2.get_p(), e2.get_valence_number(), e2.get_R0())
            calc.pp.add_pair_potential(s1,s2,vdW,eVA=False)

            # debug
            #x, y, dy = calc.pp.get_table(s1, s2)
            #import finite_differences as fd
            #dyn = fd.dx(x, y)
            #np.savetxt('pp_%s_%s.out' % ( s1, s2 ), np.transpose([x, y, dy, dyn]))

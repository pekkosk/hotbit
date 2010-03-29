#! /usr/bin/env python

import numpy as np

from ase import Atoms
from hotbit.atoms import Atoms as HotbitAtoms
from hotbit.multipole_periodicity import MultipolePeriodicity

###

NAT         = 100
SX, SY, SZ  = 10.0, 10.0, 10.0
CHARGE      = 1.0

L_MAX       = 8
N           = 5
K           = 5

###

q  = (2*np.random.random([NAT])-1)*CHARGE
q -= np.sum(q)/len(q)

a = Atoms('%iH' % NAT,
          positions  = np.random.random([NAT,3])*SX,
          charges    = q,
          cell       = [ SX, SY, SZ ],
          pbc        = True
          )

b = HotbitAtoms(
    atoms      = a,
    container  = 'Bravais'
    )

mp = MultipolePeriodicity(L_MAX, N, K)
mp.update(b, b.get_charges())



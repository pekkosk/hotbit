from ase import *
#from hotbit.wa import WedgeAtoms
#from hotbit.wa import BendAtoms
from hotbit import WedgeAtoms
from hotbit import BendAtoms

from numpy.linalg import norm as norm
#from math import pi,radians

         
# bending  test
print; print 'Bending  test... See opened windows.'
viewer = 'ag'
#viewer = 'vmd'
d = 0.5
L = 5.0
a = Atoms('H2', [[L,0,0], [L,d,0]])
a = a.repeat((3,8,3)) 
view(a, viewer=viewer)

#1
bent = BendAtoms(a, wedge_angle=pi/2)
view(bent, viewer=viewer)
m = bent.get_copies()
view(m)
m.translate((8,8,0))
m.set_cell([17,17,5])
view(m.repeat((2,3,2)))

#2a
bent = BendAtoms(a, x_bent_ceter=5.0)
#view(bent, viewer=viewer)
#m = bent.get_copies()
#view(m)
#m.translate((8,8,0))
#m.set_cell([17,17,5])
#view(m.repeat((2,3,2)))

#2b
bent2 = BendAtoms(a, x_bent_ceter=3.0)
#view(bent2, viewer=viewer)

#3
bent3 = BendAtoms(a, bend_factor=2.0)
#view(bent3, viewer=viewer)


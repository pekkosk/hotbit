from ase import *
#from hotbit.wa import WedgeAtoms
from hotbit import WedgeAtoms

from numpy.linalg import norm as norm


# get_vector test
print 'get_vector test:'
viewer = 'ag'
d = 1.5
wangle = radians(45)
t1 = radians(3)
t2 = wangle - t1 
tc = wangle / 2.0
w = WedgeAtoms(wangle, 'H3', [[d * cos(t1), d * sin(t1), 0],
                              [d * cos(tc), d * sin(tc), 0],
                              [d * cos(t2), d * sin(t2), 0]
                              ])

print 'min_d(0,2) =', norm(w.get_vector(0, 2))
print 'd(0,2)     =', norm(w.get_vector(0, 2, mic=False))

print
print 'min_d(0,1) =', norm(w.get_vector(0, 1))
print 'd(0,1)     =', norm(w.get_vector(0, 1, mic=False))

view(w, viewer=viewer)
       
# get_copies() - all    
#view(w.get_copies(), viewer=viewer)

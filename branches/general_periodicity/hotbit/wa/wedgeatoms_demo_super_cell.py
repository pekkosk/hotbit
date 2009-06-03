#
#  Super Cell Idea Demonstration
# 
from math import pi,sin,cos
from numpy import *

from ase import Atoms
from ase.visualize import *
from hotbit import WedgeAtoms

import numpy as np

# this also works!
w = WedgeAtoms(pi/2, 'OCCCCH', [[0,0,0], [2,2,0], [-2,2,0],[2,-2,0],[-2,-2,0], [0.2,0.2,0]])
w.repeat((1,1,2)).positions


d = 1
L = 2
viewer='ag'
#viewer='vmd'
#w = WedgeAtoms(pi/2, 'OCCCCH', [[0,0,0], [L,L,0], [-L,L,0],[L,-L,0],[-L,-L,0], [d,d,0]])
#w = WedgeAtoms(pi/2, 'OCH', [[0,0,0], [L,L,0], [d,d,0]])
#w = WedgeAtoms(pi/4, 'OH', [[0,0,0],[d*cos(pi/8),d*sin(pi/8),0]])
w = WedgeAtoms(radians(22.5), 'OH', [[0,0,0],[d*cos(pi/16),d*sin(pi/16),0]])
view(w, viewer=viewer)
#a = w.get_copies()
a = w.get_copies()
view(a, viewer=viewer)
a.translate([L,L,0])
a.set_cell([2*L,2*L,1])
bulk = a.repeat([3,3,5])
#bulk = a.repeat([5,5,10])
view(bulk, viewer=viewer)

#===============================================================================
# d = 1
# L = 2
# viewer='ag'
# #viewer='vmd'
# #w = WedgeAtoms(pi/2, 'OCCCCH', [[0,0,0], [L,L,0], [-L,L,0],[L,-L,0],[-L,-L,0], [d,d,0]])
# #w = WedgeAtoms(pi/2, 'OCH', [[0,0,0], [L,L,0], [d,d,0]])
# w = WedgeAtoms(pi/4, 'OH', [[0,0,0],[d*cos(pi/8),d*sin(pi/8),0]])
# #view(w, viewer=viewer)
# #a = w.get_copies(4)
# a = w.get_copies(8)
# view(a, viewer=viewer)
# 
# unitcellcenter = np.linalg.norm((a.get_cell()).forall) 
# translate_vector = a.get_center_of_mass() -  
# #for i in range(len(translate_vector)):
# #    translate_vector[i] += np.linalg.norm(a.get_cell()[i])
# a.translate(translate_vector)
# view(a, viewer=viewer)
# 
# a.set_cell([2*L,2*L,2])
# bulk = a.repeat([3,3,5])
# #bulk = a.repeat([5,5,10])
# view(bulk, viewer=viewer)
#===============================================================================


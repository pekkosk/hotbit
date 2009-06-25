from ase import *

from hotbit import WedgeAtoms
from hotbit import BendAtoms
from hotbit import GeneralizedLennardJones

from numpy.linalg import norm as norm
from numpy import Inf as numpy_Inf
#from math import pi,radians

def test(test_ok, text='Test result: '):
    if test_ok:
        print text, '\t\t-> \t\t\t Ok :)'
    else:  
        print text, '\t\t-> \t\t\t Error!'


print 'GeneralizedLennardJones test.'
viewer = 'ag'

#test structure1 (One atom)
d = 1.5
wangle1 = radians(90)
t1 = radians(0)
t2 = wangle1 - t1 
tc = wangle1 / 2.0

structure1 = [[d * cos(tc), d * sin(tc), 0]]

#test structure2 (Three atoms)
d = 3.0
wangle2 = radians(45)
t1 = radians(10)
t2 = wangle2 - t1 
tc = wangle2 / 2.0

structure2 = [[d * cos(t1), d * sin(t1), 0],
              [d / 2 * cos(tc), d / 2 * sin(tc), 0],
              [d * cos(t2), d * sin(t2), 0]
              ]

#test structure3 (Two atoms)
structure3 = [[d * cos(0), d * sin(0), 0],
              [d / 2 * cos(tc), d / 2 * sin(tc), 0]
              ]
wangle3 = wangle2



#Choose the structure
structure = structure3
wangle = wangle3

# 2nd structure control
atoms = WedgeAtoms(wangle, 'H' + str(len(structure)), structure)
atoms_ref = Atoms(atoms.get_copies(), pbc=[0, 0, 0]) 
#atoms = WedgeAtoms(wangle, 'H' + str(len(structure)), structure,
#                   number_of_cells=(2, 1, 1))
#atoms_ref = Atoms(atoms.get_copy((1)) + atoms, pbc=[0, 0, 0])


view(atoms_ref, viewer=viewer)
view(atoms, viewer=viewer)



# Compute Reference Calculator Results
#view(atoms_ref, viewer=viewer)
calc_ref = LennardJones()
atoms_ref.set_calculator(calc_ref)
#print atoms_ref.pbc
e_ref = atoms_ref.get_potential_energy()
forces_ref = atoms_ref.get_forces()


#---------------------------------------------- # Reference Calculator Results 2
#---------------------------------------------------- calc_ref2 = LennardJones()
#------------------------------------------- atoms_ref2.set_calculator(calc_ref)
#---------------------------------------------------------- #print atoms_ref.pbc
#------------------------------------ e_ref2 = atoms_ref2.get_potential_energy()
#----------------------------------------- forces_ref2 = atoms_ref2.get_forces()
#--------------------------------------------------------- print e_ref == e_ref2
#----------------------------------------------- print forces_ref == forces_ref2


# Compute New Calculator Results
#view(atoms, viewer=viewer)
#print "get_number_of_cells = ", atoms.get_number_of_cells()
calc = GeneralizedLennardJones()
atoms.set_calculator(calc)
e = atoms.get_potential_energy()
forces = atoms.get_forces()

force_sum = 0.0
for i in range(len(forces)):
    force_sum += forces[i]



# Testing, 
# 1. Energy
n_images = atoms.get_number_of_cells()[0]
#print; print 'Reference Calculator Results'
#print; print 'New Calculator Results'
print 'Energy_ref = ', e_ref
#print 'Energy_ref/', n_images, ' = ', e_ref / n_images

print 'Energy = ', e
#print 'Forces = ', forces
e2 = e_ref / n_images
#print 'Energy1 = %.20f' % e
#print 'Energy2 = %.20f' % e2
#test(e == e2)
print '|Energy1-Energy2| = %.3e' % abs(e - e2)
test(abs(e - e2) < 1e-10, 'Energy test:')
#test(abs(e - e2) < 1e-20)

# 2. Forces
print; print;
print 'Forces_ref = ', forces_ref
print 'Forces_ref/', n_images, ' = ', forces_ref[:len(structure)] / n_images
print 'Forces = ', forces
#test((forces == forces_ref/ n_images).all())
#test((abs(forces - forces_ref / n_images) < 1e-10).all(), "Forces Test")
#test((forces == forces_ref).all())

print "force_sum = ", force_sum
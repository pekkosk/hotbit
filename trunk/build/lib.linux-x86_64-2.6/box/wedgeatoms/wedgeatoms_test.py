from ase import *

from hotbit import WedgeAtoms
from hotbit import BendAtoms

from numpy.linalg import norm as norm
from numpy import Inf as numpy_Inf
#from math import pi,radians

def test(test_ok):
    if test_ok:
        print '..Ok :)'
    else:  
        print 'Error!!!'



# Comparative Atoms, WegdeAtoms test
print; print 'Comparative Atoms, WegdeAtoms test'
a = Atoms(molecule('CO'))
w = WedgeAtoms(radians(30), molecule('CO'))
print 'w.len == a.len ', w.__len__() == a.__len__()
print 'w.positions == a.positions? ', (w.positions == a.positions)
#print 'w.len = ', w.__len__()
#print 'a.len = ', a.__len__()
#print 'w.positions = ', w.positions
#print ' = ', a.positions
print 'Atoms, a:', a.__repr__()
print 'WedgeAtoms, w:', w.__repr__()
#a.pbc=True # Illegal in ASE
#print a.get_pbc() # Fails
print 'angle = ', w.angle


# __repr__ test1 (also tests set_pbc) 
print; print 'set_pbc_z, __repr__ test1 (also tests set_pbc)'
w = WedgeAtoms(radians(30), molecule('CO'))
print 'w:', w.__repr__()
w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=False)
print 'w:', w.__repr__()
w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
print 'w:', w.__repr__()

# copy() test
print; print 'copy() test. Test is ok if all are True:'
w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
new_w = w.copy()            ### Uncomment!!!!
print (new_w == w) == True   ### Uncomment!!!!
#print new_w.pbc
#print new_w.get_pbc()
#print new_w.get_pbc_z()
#print w.pbc
#w.set_pbc(False)
#print w.pbc
#print (new_w == w) == False


# __repr__() test
#print; print '__repr__() test'

# Atoms.__repr__()
#--------------------------------------- # Does not work with b = eval(repr(a))!
# # http://docs.python.org/ says: this function [repr(object)] makes an attempt to return a string that would yield an object with the same value when passed to eval(), otherwise the representation is a string enclosed in angle brackets that contains the name of the type of the object together with additional information often including the name and address of the object.
#------------------------------------------- a = Atoms(molecule('CO'), pbc=True)
#------------- b = eval(repr(a)) # results in SyntaxError due to 'positions=...'
#------------------------------------------------------------------ print v == w
#------------------------------------------ #a = Atoms(molecule('CO'), pbc=True)
#-------------------------------------------- #a.set_cell([1.0/3, 1.0/3, 1.0/3])
#---------------------------------------------------------------------- #print a

# __repr__ fails:
#w = WedgeAtoms(_pi / 6, molecule('CO'), pbc=True)
#print w
#v = eval(repr(a))
#print v
#rint v == w

#__eq__() test
print; print '__eq__() test. Test is ok if all are True:'
w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
test(w.get_pbc_z() == True)
new_w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
test((new_w == w) == True)
v = WedgeAtoms(radians(60), molecule('CO'), pbc_z=True)
u = WedgeAtoms(radians(30), molecule('CO'))
test(u.get_pbc_z() == False)
test((v == w) == False)
test((u == w) == False)
w.set_pbc_z(False)
test((new_w == w) == False) 


# get_copies, get_copy Type,ValueError tests 
viewer = 'ag'
d = 1
L = 2
#w = WedgeAtoms(radians(90/1.8), 'OCH', [[0,0,0], [L,L,0], [d,d,0]])
# view(w.get_copies(0), viewer=viewer) # ValueError
#view(w.get_copies(6), viewer=viewer) # OK!
#view(w.get_copies(7), viewer=viewer) # ValueError
#view(w.get_copies(-1,0), viewer=viewer) # ValueError
#view(w.get_copies(1,1), viewer=viewer) # ValueError
#view(w.get_copies(1,12), viewer=viewer) # ValueError
#view(w.get_copy((3,4)), viewer=viewer) # TypeError

# get_copies, get_copy test
viewer = 'ag'
d = 1
L = 2
w = WedgeAtoms(radians(90), 'OCH', [[0, 0, 0], [L, L, 0], [d, d, 0]])
#view(w.get_copies(), viewer=viewer)
#view(w.get_copies(1,-3), viewer=viewer)
#view(w.get_copies(2), viewer=viewer)
#view(w.get_copies(0,-2), viewer=viewer)
#view(w.get_copy(0.5), viewer=viewer)

#w = WedgeAtoms(pi/6, 'OC', [[0,0,0], [0.3,0.2,0.]])
#w = WedgeAtoms(pi/4, 'OC', [[0,0,0], [0.3,0.3,0.]])
#w = WedgeAtoms(pi/4, 'OCCCCH', [[0,0,0], [0.3,0.3,0.], [0.3,-0.3,0.],[-0.3,0.3,0.],[-0.3,-0.3,0.], [0.2,0.1,0]])
#view(w.get_copies(8), viewer=viewer)


# get_vector and get_distance test
print; print 'get_vector test:'
viewer = 'ag'
d = 1.5
wangle = radians(45)
t1 = radians(3)
t2 = wangle - t1 
tc = wangle / 2.0
w = WedgeAtoms(wangle, 'H3', [[d * cos(t1), d * sin(t1), 0],
                              [d * cos(tc), d * sin(tc), 0],
                              [d * cos(t2), d * sin(t2), 0]
                              ], cell=(0, 0, 10))
print 'min_d(0,2) =', norm(w.get_vector(0, 2))
print 'd(0,2)     =', norm(w.get_vector(0, 2, radial_mic=False))
print
print 'min_d(0,1) =', norm(w.get_vector(0, 1))
print 'd(0,1)     =', norm(w.get_vector(0, 1, radial_mic=False))
#view(w, viewer=viewer)

# get_distance test
print; print 'get_distance test:'
print
for i in range(len(w)):
    for j in range(len(w)):   
        if w.get_distance(i, j, mic=True) != w.get_distance(i, j, mic=False):
            print 'norm of get_vector(%d,%d) MIC =' % (i, j), w.get_distance(i, j, mic=True)
            print 'norm of get_vector(%d,%d)     =' % (i, j), w.get_distance(i, j, mic=False)
            



# See internal pbc representation
print; print 'See internal pbc representation'
w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
#print 'WedgeAtoms, w:', w.__repr__()
print 'w.get_pbc() = ', w.get_pbc()
w.set_pbc(False)
print 'w.get_pbc() = ', w.get_pbc()
w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
print 'w.pbc_z = ', w.pbc_z

 
# set_pbc_z, get_pbc_z test 
print; print 'set_pbc_z, get_pbc_z test.'
w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
test(w.get_pbc_z() == True)
w.set_pbc_z(False)
test(w.get_pbc_z() == False)
w.set_pbc_z(True)
test(w.get_pbc_z() == True)

#1
a = Atoms(molecule('CO'))
w = WedgeAtoms(radians(30), a)
test(w.get_pbc_z() == False)
#print w

#2
a = Atoms(molecule('CO'), pbc=True)
w = WedgeAtoms(radians(30), a)
test(w.get_pbc_z() == True)  
#print w


# set_pbc of bad inputs, generates TypeError
#w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
#w.set_pbc_z([False,False,False]) # TypeError
#w.set_pbc((False,False,False)) # TypeError

#print; print "Example of Illegal operations"
#w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
#print 'Is pbc accessible?', w.pbc   # NOT Fails!
#w.pbc = False # Illegal in ASE
#print w.get_pbc() 
#w.pbc = True
#print w.get_pbc() # Fails


     

#Atoms -> WedgeAtoms
print; print "#Atoms -> WedgeAtoms"
w = WedgeAtoms(radians(30))
a = Atoms(w)
# FIXME: this might be not ok!
print "Default pbc then are: ", a.get_pbc() , '-- surprise, surprise!'




# cell test
print; print 'cell test.. see opened windows, if any'
viewer = 'ag'
d = 1.5
wangle = radians(45)
t1 = radians(3)
t2 = wangle - t1 
tc = wangle / 2.0
structure = [[d * cos(t1), d * sin(t1), 0],
             [d * cos(tc), d * sin(tc), 0],
             [d * cos(t2), d * sin(t2), 0]
             ]

w = WedgeAtoms(wangle, 'H3', structure)
#view(w, viewer=viewer)

w = WedgeAtoms(wangle, 'H3', structure, cell=10)
#view(w, viewer=viewer)

w = WedgeAtoms(wangle, 'H3', structure, cell=(10,2))
#view(w, viewer=viewer)

w = WedgeAtoms(wangle, 'H3', structure, cell=(2, 10))
#view(w, viewer=viewer)

w = WedgeAtoms(wangle, 'H3', structure, cell=(2, (0, 0, 10)))
#view(w, viewer=viewer)

#NonImplementedError
#w = WedgeAtoms(wangle, 'H3', structure, cell=(2, (0, 10, 10)))  
#view(w, viewer=viewer)

#w = WedgeAtoms(wangle, 'H3', structure, cell=('er',)) #ValueError
#view(w, viewer=viewer)

#w = WedgeAtoms(wangle, 'H3', structure, cell=(2, 'er')) #TypeError
#view(w, viewer=viewer)

#w = WedgeAtoms(wangle, 'H3', structure, cell='er') #TypeError
#view(w, viewer=viewer)



# cell daring test
print; print 'cell daring test.. see opened windows, if any'
viewer = 'ag'
w = WedgeAtoms(wangle, 'H3', structure, cell=3)
#view(w, viewer=viewer)



# symmetry_operation test
print; print "symmetry_operation test:"
d = 1
deta_z = 5
structure = [[0, 0, 0], [d * cos(pi / 16), d * sin(pi / 16), 0]]
w = WedgeAtoms(radians(22.5), 'OH', structure, cell=(d,deta_z))
#view(w) 
#w.symmetry_operation(1,(20,3)) # ValueError
w.set_number_of_cells((5, numpy_Inf))
#w.symmetry_operation(1,(6,3)) # ValueError
#w.symmetry_operation(1,(-2,3)) # ValueError
#w.symmetry_operation(1,(2,3,1)) # ValueError 
#w.symmetry_operation(1,10) # TypeError
print '1'
print w.symmetry_operation(1, (0, 3, 0))
print w.get_vector(1,0)
print w.symmetry_operation(1, (0, 0))
print w.symmetry_operation(1, (0, 2))
print w.symmetry_operation(1, (2, 0))
print w.symmetry_operation(1, (2, 2))
print w.get_copy(2).positions
print w.symmetry_operation(1, (2, -2))
# FIXME: properly check this!

#print '2'
#print w.symmetry_operation2(1, (0, 3, 0))
#print w.get_vector(1,0)
#print w.symmetry_operation2(1, (0, 0))
#print w.symmetry_operation2(1, (0, 2))
#print w.symmetry_operation2(1, (2, 0))
#print w.symmetry_operation2(1, (2, 2))
#print w.get_copy(2).positions
#print w.symmetry_operation2(1, (2, -2))
## FIXME: properly check this!



d = 1.5
wangle = radians(45)
t1 = radians(3)
t2 = wangle - t1 
tc = wangle / 2.0
structure = [[d * cos(t1), d * sin(t1), 0],
             [d * cos(tc), d * sin(tc), 0],
             [d * cos(t2), d * sin(t2), 0]
             ]

w = WedgeAtoms(wangle, 'H3', structure, cell=(30,2))
#view(w)
#view(w.get_copies(2,0))
#view(w.get_copies(2,0).repeat((1,1,4)))
#FIXME: repeat!!
#view(w.get_copies(2,0).repeat((3,2,4)))
for index in range(len(w)):
    print 'atom index:', index
    #print w.symmetry_operation(index, (0, 0))
    #print w.symmetry_operation(index, (1, 0))
    #print w.symmetry_operation(index, (0, 1))
    #print w.symmetry_operation(index, (1, 1))
    print w.symmetry_operation(index, (0, 0))
    print w.symmetry_operation(index, (2, 0))
    print w.symmetry_operation(index, (0, -2))
    print w.symmetry_operation(index, (2, -2))
    
    


# repeat test
print; print "repeat test:"

d = 1.5
wangle = radians(45)
t1 = radians(3)
t2 = wangle - t1 
tc = wangle / 2.0
structure = [[d * cos(t1), d * sin(t1), 0],
             [d * cos(tc), d * sin(tc), 0],
             [d * cos(t2), d * sin(t2), 0]
             ]

w = WedgeAtoms(wangle, 'H3', structure, cell=(30,2))
#view(w.get_copies(2,0).repeat((1,1,4)))
# FIXME: rewrite repeat and copy!!!! 
#view(w.repeat((2,2,1)))







# Bending  test
print; print 'bending  test... See opened windows, if any'
viewer = 'ag'
#viewer = 'vmd'
d = 0.5
L = 5.0
a = Atoms('H2', [[L, 0, 0], [L, d, 0]])
a = a.repeat((3, 8, 3)) 
#L = 6.0
#view(a, viewer=viewer)

#bent = BendAtoms(a, x_bent_ceter=5.0)
#view(bent, viewer=viewer)
#m = bent.get_copies()
#view(m)
#m.translate((8,8,0))
#m.set_cell([17,17,5])
#view(m.repeat((2,3,2)))

#bent2 = BendAtoms(a, x_bent_ceter=3.0)
#view(bent2, viewer=viewer)


print "# Before::: bent = BendAtoms(a, wedge_angle=pi/2)"
bent = BendAtoms(a, wedge_angle=pi / 2)
#view(bent, viewer=viewer)
m = bent.get_copies()
#view(m)
#m.translate((8,8,0))
#m.set_cell([17,17,5])
#view(m.repeat((2,3,2)))

#bent3 = BendAtoms(a, bend_factor=2.0)
#view(bent3, viewer=viewer)



# more wired banding...
d = 0.5
L = 5.0
a = Atoms('H3', [[-L, 0, 0], [-L, d, 0], [-L, -d, 0]])
#view(a, viewer=viewer)
#a +=
#b = a.copy()
#b = a + a.translate((0,0,5))  
#view(b, viewer=viewer)
a = a.repeat((3, 8, 3))
#view(a, viewer=viewer)
bent = BendAtoms(a, wedge_angle=pi / 2)
#view(bent, viewer=viewer)
#m = bent.get_copies()

#d = 0.5
#L = 5.0
#a = Atoms('H2', [[L, 0, 0], [L, -d, 0]])
#a = a.repeat((3, 8, 3))
#view(a, viewer=viewer)
#bent = BendAtoms(a, wedge_angle=pi / 2)
#view(bent, viewer=viewer)




# arrays test
for name, a in w.arrays.items():
    print name, ' = ', a
    
# number_of_cells test
print; print "number_of_cells test:"
w = WedgeAtoms(radians(90), pbc_z=True)
print w.get_number_of_cells()
w = WedgeAtoms(radians(90))
print w.get_number_of_cells()
#w.set_number_of_cells((-1,2,2)) # ValueError
#w.set_number_of_cells((1,-2,2)) # ValueError
#w.set_number_of_cells((1,2,-3)) # ValueError
#w.set_number_of_cells((1,2,3,4)) # ValueError
# w.set_number_of_cells((1,1,2)) # ValueError
#w.set_number_of_cells((5,1)) # ValueError
w.set_number_of_cells((2, 1))
w.set_number_of_cells((2, numpy_Inf))
print w.get_number_of_cells()
w = WedgeAtoms(radians(90), pbc_z=True,
               number_of_cells=(2, numpy_Inf))
print w



# transfoms_der test
print; print "transfoms_der test"
print "Rotations for 45 degrees"
w = WedgeAtoms(radians(45), molecule('CO')) 
print w.transform_der((0,0,0))
print w.transform_der((1,0,0))
print w.transform_der((-1,3,0))
print w.transform_der((2,0,0))
print w.transform_der((3,0,0))
print w.transform_der((4,0,0))
print w.transform_der((8,0,0))
print "Rotations for n*30 degrees"
w = WedgeAtoms(radians(30), molecule('CO')) 
print w.transform_der((0,0,0))
print w.transform_der((1,0,0))
print w.transform_der((2,0,0))
print w.transform_der((3,0,0))
from ase import *
#from hotbit.wa import WedgeAtoms
#from hotbit.wa import BendAtoms
from hotbit import WedgeAtoms
from hotbit import BendAtoms

from numpy.linalg import norm as norm
#from math import pi,radians

# Pure set_pbc test: See internal pbc representation
print; print 'Pure set_pbc test: See internal pbc representation'
w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
#print 'WedgeAtoms, w:', w.__repr__()
print w.get_pbc()
w.set_pbc(False)
print w.get_pbc()
w.set_pbc(True)
print w.get_pbc()


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
new_w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
print (new_w == w) == True
v = WedgeAtoms(radians(60), molecule('CO'), pbc_z=True)
u = WedgeAtoms(radians(30), molecule('CO'))
print (v == w) == False
print (u == w) == False
w.set_pbc_z(False)
print (new_w == w) == False

# arrays test
for name, a in w.arrays.items():
    print name, ' = ', a 

# set_pbc_z, get_pbc_z test 
print; print 'set_pbc_z, get_pbc_z test. Test is ok if all are True:'
w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
print w.get_pbc_z() == True
w.set_pbc_z(False)
print w.get_pbc_z() == False
w.set_pbc_z(True)
print w.get_pbc_z() == True

# set_pbc_z, __repr__ test1 (also tests set_pbc) 
print; print 'set_pbc_z, __repr__ test1 (also tests set_pbc)'
w = WedgeAtoms(radians(30), molecule('CO'))
print 'w:', w.__repr__()
w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=False)
print 'w:', w.__repr__()
w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
print 'w:', w.__repr__()


w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
print w.pbc_z

print; print 'Is pbc accessible?', w.pbc   # NOT Fails!

#Atoms -> WedgeAtoms
a = Atoms(w)
print a, 'Hohoohhhhhhh!'
 

# set_pbc of bad inputs
#print; print 'set_pbc test2'
#w = WedgeAtoms(radians(30), molecule('CO'), pbc_z=True)
#w.set_pbc_z([False,False,False])
#w.set_pbc((False,False,False))
#print w.get_pbc()

#w.pbc = False # Illegal in ASE
#print w.get_pbc() 
#w.pbc = True
#rint w.get_pbc() # Fails


# this also works!
w.repeat((2, 1, 2)).positions

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
                              ])
print 'min_d(0,2) =', norm(w.get_vector(0, 2))
print 'd(0,2)     =', norm(w.get_vector(0, 2, mic=False))
print
print 'min_d(0,1) =', norm(w.get_vector(0, 1))
print 'd(0,1)     =', norm(w.get_vector(0, 1, mic=False))
#view(w, viewer=viewer)

# get_distance test
print; print 'get_distance test:'
print
for i in range(len(w)):
    for j in range(len(w)):   
        if w.get_distance(i, j, mic=True) != w.get_distance(i, j, mic=False):
            print 'norm of get_vector(%d,%d) MIC =' % (i, j), w.get_distance(i, j, mic=True)
            print 'norm of get_vector(%d,%d)     =' % (i, j), w.get_distance(i, j, mic=False)
            

# Bending  test
print; print 'bending  test... See opened windows.'
viewer = 'ag'
#viewer = 'vmd'
d = 0.5
L = 5.0
a = Atoms('H2', [[L,0,0], [L,d,0]])
a = a.repeat((3,8,3)) 
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

bent = BendAtoms(a, wedge_angle=pi/2)
view(bent, viewer=viewer)
#m = bent.get_copies()
#view(m)
#m.translate((8,8,0))
#m.set_cell([17,17,5])
#view(m.repeat((2,3,2)))

#bent3 = BendAtoms(a, bend_factor=2.0)
#view(bent3, viewer=viewer)


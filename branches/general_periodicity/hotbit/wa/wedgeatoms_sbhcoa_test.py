from math import pi as _pi

from ase import *
#from hotbit.wa import WedgeAtoms
from hotbit.wa.wedgeatoms_sbhcoa import WedgeAtoms

# Pure set_pbc test: See internal pbc representation
print; print 'Pure set_pbc test: See internal pbc representation'
w = WedgeAtoms(_pi / 6, molecule('CO'), pbc_z=True)
#print 'WedgeAtoms, w:', w.__repr__()
print w.get_pbc()
w.set_pbc(False)
print w.get_pbc()
w.set_pbc(True)
print w.get_pbc()


# Comparative Atoms, WegdeAtoms test
print; print 'Comparative Atoms, WegdeAtoms test'
a = Atoms(molecule('CO'))
w = WedgeAtoms(_pi / 6, molecule('CO'))
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
print 'angle = ', w.get_angle()



# copy() test
print; print 'copy() test. Test is ok if all are True:'
w = WedgeAtoms(_pi / 6, molecule('CO'), pbc_z=True)
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
w = WedgeAtoms(_pi / 6, molecule('CO'), pbc_z=True)
new_w = WedgeAtoms(_pi / 6, molecule('CO'), pbc_z=True)
print (new_w == w) == True
v = WedgeAtoms(_pi / 3, molecule('CO'), pbc_z=True)
u = WedgeAtoms(_pi / 6, molecule('CO'))
print (v == w) == False
print (u == w) == False
w.set_pbc_z(False)
print (new_w == w) == False

# arrays test
for name, a in w.arrays.items():
    print name, ' = ', a 

# set_pbc_z, get_pbc_z test 
print; print 'set_pbc_z, get_pbc_z test. Test is ok if all are True:'
w = WedgeAtoms(_pi/6, molecule('CO'), pbc_z=True)
print w.get_pbc_z() == True
w.set_pbc_z(False)
print w.get_pbc_z() == False
w.set_pbc_z(True)
print w.get_pbc_z() == True

# set_pbc_z, __repr__ test1 (also tests set_pbc) 
print; print 'set_pbc_z, __repr__ test1 (also tests set_pbc)'
w = WedgeAtoms(_pi / 6, molecule('CO'))
print 'w:', w.__repr__()
w = WedgeAtoms(_pi / 6, molecule('CO'), pbc_z=False)
print 'w:', w.__repr__()
w = WedgeAtoms(_pi / 6, molecule('CO'), pbc_z=True)
print 'w:', w.__repr__()


w = WedgeAtoms(_pi / 6, molecule('CO'), pbc_z=True)
print w.pbc_z
#print w.pbc         # Fails!!!!! =)


# set_pbc of bad inputs
#print; print 'set_pbc test2'
#w = WedgeAtoms(_pi / 6, molecule('CO'), pbc_z=True)
#w.set_pbc_z([False,False,False])
#w.set_pbc((False,False,False))
#print w.get_pbc()

#w.pbc = False # Illegal in ASE
#print w.get_pbc() 
#w.pbc = True
#rint w.get_pbc() # Fails
import numpy as nu
from ceig import ceigr, ceigc
from hotbit.fortran.eigensolver import geig, geigc

A=nu.array(((1,4),(4,7)), nu.float)
B=nu.array(((1,0.2),(0.2,1)), nu.float)

a = A.copy()
b = B.copy()

ev, vec = geig(A,B)
ev_, vec_ = ceigr(a,b)

print vec
print vec_
assert nu.sum(abs(ev - ev_)) < 1e-6
assert nu.sum(abs(vec - vec_)) < 1e-6
print "Real OK"

A = nu.array(((1, 4+0.1j),(4-0.1j, 7)), nu.complex, order='c')
B = nu.array(((1,0.2+0.2j),(0.2-0.2j,1)), nu.complex)
a = A.copy()
b = B.copy()

ev, vec = geigc(A,B)
ev_, vec_ = ceigc(a,b)

print vec
print vec_
assert nu.sum(abs(ev - ev_)) < 1e-6
assert nu.sum(abs(vec - vec_)) < 1e-6
print "Complex OK"

import numpy as nu
from ceig import ceigr, ceigc
from hotbit.fortran.eigensolver import geig, geigc

A=nu.array(((1,4),(4,7)), nu.float)
B=nu.array(((1,0.2),(0.2,1)), nu.float)

a = A.copy()
b = B.copy()

ev, vec = geig(A,B)
ev_, vec_ = ceigr(a,b)

#print vec
#print vec_
err_ev = nu.sum(abs(ev - ev_))
err_vec = nu.sum(abs(vec - vec_))
if err_ev < 1e-15 and err_vec < 1e-15:
    print "Real generalized eigensolver OK"
    print "    error(eigenvalues): %f" % err_ev
    print "    error(eigenvectors): %f" % err_vec
    print ""

A = nu.array(((1, 4+0.1j),(4-0.1j, 7)), nu.complex, order='c')
B = nu.array(((1,0.2+0.2j),(0.2-0.2j,1)), nu.complex)
a = A.copy()
b = B.copy()

ev, vec = geigc(A,B)
ev_, vec_ = ceigc(a,b)

#print vec
#print vec_
err_ev = nu.sum(abs(ev - ev_))
err_vec = nu.sum(abs(vec - vec_))
if err_ev < 1e-15 and err_vec < 1e-15:
    print "Complex generalized eigensolver OK"
    print "    error(eigenvalues): %f" % err_ev
    print "    error(eigenvectors): %f" % err_vec
    print ""
#assert nu.sum(abs(ev - ev_)) < 1e-15
#assert nu.sum(abs(vec - vec_)) < 1e-15
#print "Complex OK"

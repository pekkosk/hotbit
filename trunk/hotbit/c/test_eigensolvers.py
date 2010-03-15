import numpy as nu
from _hotbit import geig as geig_c
from hotbit.fortran.eigensolver import geig as geig_for
from hotbit.fortran.eigensolver import geigc as geigc_for

THRES=1e-14

A=nu.array(((1,4),(4,7)), nu.float)
B=nu.array(((1,0.2),(0.2,1)), nu.float)

a = A.copy()
b = B.copy()

ev, vec = geig_for(A,B)
ev_, vec_ = geig_c(a,b)

#print vec
#print vec_
err_ev = nu.sum(abs(ev - ev_))
err_vec = nu.sum(abs(vec - vec_))
if err_ev < THRES and err_vec < THRES:
    print "Real generalized eigensolver OK"
else:
    print "Real generalized eigensolver FAILED"

print "    error(eigenvalues): %e" % err_ev
print "    error(eigenvectors): %e" % err_vec
print ""
    
A = nu.array(((1, 4+0.1j),(4-0.1j, 7)), nu.complex, order='c')
B = nu.array(((1,0.2+0.2j),(0.2-0.2j,1)), nu.complex)
a = A.copy()
b = B.copy()

ev, vec = geigc_for(A,B)
ev_, vec_ = geig_c(a,b)

#print vec
#print vec_
err_ev = nu.sum(abs(ev - ev_))
err_vec = nu.sum(abs(vec - vec_))
if err_ev < THRES and err_vec < THRES:
    print "Complex generalized eigensolver OK"
else:
    print "Complex generalized eigensolver FAILED"

print "    error(eigenvalues): %e" % err_ev
print "    error(eigenvectors): %e" % err_vec
print ""
#assert nu.sum(abs(ev - ev_)) < 1e-15
#assert nu.sum(abs(vec - vec_)) < 1e-15
#print "Complex OK"

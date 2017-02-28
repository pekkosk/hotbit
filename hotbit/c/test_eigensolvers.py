import numpy as np
from _hotbit import geig as geig_c
from hotbit.c.old_fortran.eigensolver import geig as geig_for
from hotbit.c.old_fortran.eigensolver import geigc as geigc_for

THRES=1e-14

A=np.array(((1,4),(4,7)), np.float)
B=np.array(((1,0.2),(0.2,1)), np.float)

a = A.copy()
b = B.copy()

ev, vec = geig_for(A,B)
ev_, vec_ = geig_c(a,b)

#print vec
#print vec_
err_ev = np.sum(abs(ev - ev_))
err_vec = np.sum(abs(vec - vec_))
if err_ev < THRES and err_vec < THRES:
    print("Real generalized eigensolver OK")
else:
    print("Real generalized eigensolver FAILED")

print("    error(eigenvalues): %e" % err_ev)
print("    error(eigenvectors): %e" % err_vec)
print("")
    
A = np.array(((1, 4+0.1j),(4-0.1j, 7)), np.complex, order='c')
B = np.array(((1,0.2+0.2j),(0.2-0.2j,1)), np.complex)
a = A.copy()
b = B.copy()

ev, vec = geigc_for(A,B)
ev_, vec_ = geig_c(a,b)

#print vec
#print vec_
err_ev = np.sum(abs(ev - ev_))
err_vec = np.sum(abs(vec - vec_))
if err_ev < THRES and err_vec < THRES:
    print("Complex generalized eigensolver OK")
else:
    print("Complex generalized eigensolver FAILED")

print("    error(eigenvalues): %e" % err_ev)
print("    error(eigenvectors): %e" % err_vec)
print("")
#assert np.sum(abs(ev - ev_)) < 1e-15
#assert np.sum(abs(vec - vec_)) < 1e-15
#print "Complex OK"

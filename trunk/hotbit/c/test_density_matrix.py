from hotbit.fortran.misc import fortran_rhoc
from density_matrix import c_rhoc
import numpy as nu

norb = 3
nk = 1
wf = nu.random.random((nk,norb,norb)) + 1j*nu.random.random((nk,norb,norb))
occ = nu.zeros((nk, norb), nu.float)
occ[:,0:norb/2+1] = 2.0

rho_f = fortran_rhoc(wf, occ, norb, nk)
rho_c = c_rhoc(wf, occ)

print rho_f
print rho_c
if nu.sum(abs(rho_f-rho_c)) < 1e-15:
    print "Complex density matrix OK"

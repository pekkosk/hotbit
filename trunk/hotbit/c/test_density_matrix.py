from hotbit.fortran.misc import fortran_rhoc, fortran_rhoec
from density_matrix import c_rhoc, c_rhoec
import numpy as nu

norb = 5
nk = 3
wf = nu.random.random((nk,norb,norb)) + 1j*nu.random.random((nk,norb,norb))
occ = nu.zeros((nk, norb), nu.float)
occ[:,0:norb/2+1] = 2.0
e = nu.random.random((nk, norb))

rho_f = fortran_rhoc(wf, occ, norb, nk)
rho_c = c_rhoc(wf, occ)

#print rho_f
#print rho_c
err = nu.sum(abs(rho_f-rho_c))
if err < 1e-15:
    print "Complex density matrix OK"
    print "    error: %f" % err
    print ""

rhoe_f = fortran_rhoec(wf, occ, e, norb, nk)
rhoe_c = c_rhoec(wf, occ, e)
#print rhoe_f
#print rhoe_c
err = nu.sum(abs(rhoe_f-rhoe_c))
if err < 1e-15:
    print "Energy-weighted complex density matrix OK"
    print "    error: %f" % err
    print ""



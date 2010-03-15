from hotbit.fortran.misc import fortran_rhoc, fortran_rhoec
from _hotbit import c_rhoc, c_rhoec
import numpy as nu

THRES = 1e-12

def py_rho(wf, occ):
    nk, n, m = wf.shape
    nk2, k = occ.shape

    assert n == m
    assert n == k
    assert nk == nk2

    rho = nu.zeros(wf.shape, dtype=wf.dtype)
    for k in range(nk):
        rho[k,:,:] = nu.dot(
            wf[k,:,:].transpose(),
            occ[k,:].reshape(n,-1) * wf[k,:,:].conj()
            )
    return rho


def py_rhoe(wf, occ, e):
    nk, n, m = wf.shape
    nk2, k = occ.shape

    assert n == m
    assert n == k
    assert nk == nk2

    rhoe = nu.zeros(wf.shape, dtype=wf.dtype)
    for k in range(nk):
        rhoe[k,:,:] = nu.dot(
            wf[k,:,:].transpose(),
            e[k,:].reshape(n,-1)*occ[k,:].reshape(n,-1) * wf[k,:,:].conj()
            )
    return rhoe
    

norb = 5
nk = 3
wf = nu.random.random((nk,norb,norb)) + 1j*nu.random.random((nk,norb,norb))
#for k in range(nk):
#    wf[k,:,:] = (wf[k,:,:]+nu.transpose(wf[k,:,:]))/2
occ = nu.zeros((nk, norb), nu.float)
occ[:,0:norb/2+1] = 2.0
e = nu.random.random((nk, norb))

rho_f = fortran_rhoc(wf, occ, norb, nk)
rho_c = c_rhoc(wf, occ)
rho_p = py_rho(wf, occ)

err = nu.sum(abs(rho_f-rho_c))
if err < 1e-15:
    print "Complex density matrix OK"
print "    error: %e" % err
print ""

err = nu.sum(abs(rho_f-rho_p))
if err < THRES:
    print "Complex density matrix OK (Python)"
print "    error: %e" % err
print ""

rhoe_f = fortran_rhoec(wf, occ, e, norb, nk)
rhoe_c = c_rhoec(wf, occ, e)
rhoe_p = py_rhoe(wf, occ, e)
#print rhoe_f
#print rhoe_c
err = nu.sum(abs(rhoe_f-rhoe_c))
if err < THRES:
    print "Energy-weighted complex density matrix OK"
print "    error: %e" % err
print ""

err = nu.sum(abs(rhoe_f-rhoe_p))
if err < THRES:
    print "Energy-weighted complex density matrix OK (Python)"
print "    error: %e" % err
print ""




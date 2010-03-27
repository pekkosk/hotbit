from ase import *
from pylab import *
from ase import Atoms as ase_Atoms
from hotbit import *
from hotbit.atoms import Atoms
from hotbit.test.misc import default_param
from numpy import *

eps = 1E-6


# check that C1H1-presentation of C6H6 goes right
n1 = 4
SCC=True
cut=1E10
a1 = Atoms('CH',[(1.42,0,0),(2.42,0.1,0.1)],container='Wedge')
a1.set_container(M=6,height=10)
aux = a1.extended_copy([(1,0,0)])
a1 += aux
a1.set_container(M=3,height=10)


n2 = 12
c1 = Hotbit(SCC=SCC,txt='-',kpts=(3,1,1),gamma_cut=cut,**default_param)
a1.set_calculator(c1)
a1.get_potential_energy()


a2 = a1.extended_copy((3,1,1))
c2 = Hotbit(SCC=SCC,txt='-',kpts=(1,1,1),gamma_cut=cut,**default_param)
a2.set_calculator(c2)
a2.get_potential_energy()

# start Mulliken analysis
q1 = c1.get_dq()
q2 = c2.get_dq()[:4]
assert all(abs(q1-q2)<eps)

# atom mulliken
for i in range(4):
    q1 = c1.get_atom_mulliken(i)
    q2 = c2.get_atom_mulliken(i)
    assert abs(q1-q2)<eps
    
# basis mulliken
for mu in [0,4,5,9]: #s-orbital occupations should be the same
    q1 = c1.get_basis_mulliken(mu)
    q2 = c2.get_basis_mulliken(mu)
    assert abs(q1-q2)<eps
    
# atom wf mulliken
norb = c1.st.norb
q = zeros(4)
for i in range(4):
    q1 = 0.0
    for a in range(norb):
        for k in range(c1.st.nk):        
            occ = c1.st.f[k,a]
            wk = c1.st.wk[k]
            q1 += c1.get_atom_wf_mulliken(i,k,a,wk=False) * occ * wk
    q2 = c1.get_atom_mulliken(i)
    assert abs(q2-q1)<eps

norb = c2.st.norb
for i in range(12):
    q1 = 0.0
    for a in range(norb):
        q1 += c2.get_atom_wf_mulliken(i,0,a)*c2.st.f[0,a]
    q2 = c2.get_atom_mulliken(i)
    assert abs(q2-q1)<eps
    
    
# atom all angmom
i=0
q = zeros((3,))
for a in range(c2.st.norb):
    q += c2.get_atom_wf_all_angmom_mulliken(i,0,a)
assert all((q-(1,3,0))<eps)

i=0
q = zeros((3,))
for a in range(c1.st.norb):
    for k in range(c1.st.nk):
        q += c1.get_atom_wf_all_angmom_mulliken(i,k,a)
assert all((q-(1,3,0))<eps)


# total DOS
if False:
    e1,dos1 = c1.get_density_of_states(broaden=True)
    e2,dos2 = c2.get_density_of_states(broaden=True)
    plot(e1,dos1,label='4 atoms')
    plot(e2,dos2/3,label='12 atoms')
    legend()
    show()


npts = 2000
e2,ldos2 = c2.get_local_density_of_states(width=0.1,npts=npts,window=(-20,35))
e1,ldos1 = c1.get_local_density_of_states(width=0.1,npts=npts,window=(-20,35))

if False:
    for i in range(4):
        plot(e1,ldos1[i]+0.5*i)
        plot(e2,ldos2[i]+0.5*i)
    show()
    
# for all atoms the integrated LDOS should equal the number of orbitals
grid = len(e1)
de = e1[1]-e1[0]
q = zeros((4))
for i in range(4):
    for g in range(grid):
        q[i] += ldos1[i,g]*de
assert all( abs(q-(4,1,4,1))<eps )
    
    

# PLDOS (sum over angular momenta should give LDOS)
e,ldos = c1.get_local_density_of_states()
e,ldos_, pldos = c1.get_local_density_of_states(projected=True)
ldos2 = pldos.sum(axis=1)
assert all( abs(ldos-ldos2)<eps )



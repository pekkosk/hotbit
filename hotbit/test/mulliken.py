from ase import *
try:
    from pylab import *
except:
    pass
from ase import Atoms as ase_Atoms
from hotbit import *
from hotbit.atoms import Atoms
from hotbit.test.misc import default_param
from numpy import *
from box.systems import graphene
from ase.build import molecule

eps = 1E-5

if True:
    atoms = graphene(2,2,1.42)
    calc = Hotbit(SCC=False,kpts=(4,4,1),txt='-',**default_param)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    
    # compare to Fig.4 of Koskinen & Makinen, CMS 47, 237 (2009)
    if False:
        w=0.2
        x,y = calc.get_covalent_energy('default',width=w,window=(-20,15))
        plot(x,y,label='total')
        axhline(0)
        x,y = calc.get_covalent_energy('angmom',i=0,j=0,width=w,window=(-20,15))
        plot(x,y,label='ss')
        x,y = calc.get_covalent_energy('angmom',i=0,j=1,width=w,window=(-20,15))
        plot(x,y,label='sp')
        x,y = calc.get_covalent_energy('angmom',i=1,j=1,width=w,window=(-20,15))
        plot(x,y,label='pp')
        legend()
        show()
    
    
    x0,y0 = calc.get_covalent_energy('default')
    x1,y1 = calc.get_covalent_energy('angmom',i=0,j=0)
    x2,y2 = calc.get_covalent_energy('angmom',i=0,j=1)
    x3,y3 = calc.get_covalent_energy('angmom',i=1,j=1)
    assert abs(y0.sum()-(y1+y2+y3).sum())<eps
    
    atoms = graphene(1,1,1.42)
    calc = Hotbit(SCC=False,kpts=(1,1,1),txt='-',**default_param)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    
    x,y = calc.get_covalent_energy('default')
    x1,y1 = calc.get_covalent_energy('atoms',i=0,j=0)
    x2,y2 = calc.get_covalent_energy('atoms',i=0,j=1)
    x3,y3 = calc.get_covalent_energy('atoms',i=1,j=1)
    assert abs((y1+y2+y3).sum()-y.sum())<eps
    
    sm = 0.0
    for i in range(calc.st.norb):
        for j in range(i,calc.st.norb):
            x1,y1 = calc.get_covalent_energy('orbitals',i=i,j=j) 
            sm += y1.sum()
    assert abs(y.sum()-sm)<eps


if True:
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
    
    





#
#   BONDING ANALYSIS
#

#
#   Mayer bond order
#
if True:
    atoms = graphene(4,4,1.42)
    calc = Hotbit(SCC=False,kpts=(8,8,1),txt='-',**default_param)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    assert abs(calc.get_mayer_bond_order(1,2)-1.24155188722)<eps
    assert abs(calc.get_atom_energy(0)-6.95260830265)<eps
    assert abs(calc.get_atom_and_bond_energy(0)--9.62628777865)<eps
    assert abs(calc.get_promotion_energy(0)-6.95260830265)<eps
    assert abs(calc.get_bond_energy(1,2)--12.0027553172)<eps

    #print 'graphene'
    #print 'A_C:',calc.get_atom_energy(0)
    #print 'AB_C:',calc.get_atom_and_bond_energy(0)
    #print 'prom_C',calc.get_promotion_energy(0)
    #print 'B_CC:',calc.get_bond_energy(1,2)


if True:
    #
    # Analyzing absolute atom and bond energies
    # 1) using CH-modeled benzene
    # 2) using full C6H6
    #
    ## using CH to model benzene
    a = Atoms('CH',[(1.395,0,0),(2.482,0,0)],container='Wedge')
    a.container.set(M=6,height=3)
    calc = Hotbit(SCC=False,kpts=(6,1,1),txt='-',**default_param)
    a.set_calculator(calc)
    e1 = a.get_potential_energy()
    
    promC = calc.get_promotion_energy(0)
    promH = calc.get_promotion_energy(1)
    AC = calc.get_atom_energy(0)
    AH = calc.get_atom_energy(1)
    # bond energy exists here ONLY between C and H
    BCH = calc.get_bond_energy(0,1)
    ABC = calc.get_atom_and_bond_energy(0)
    ABH = calc.get_atom_and_bond_energy(1)
    
    #print 'CH'
    #print 'prom',promC,promH
    #print 'atom',AC, AH
    #print 'bond',BCH
    #print 'ab',ABC,ABH
    assert abs(e1-(AC+AH+BCH))<eps
    assert abs(e1-sum([calc.get_atom_and_bond_energy(i) for i in range(2)] ))<eps
       
    # 2) Using the full benzene
    atoms = a.extended_copy(((-2,3),1,1))
    calc = Hotbit(SCC=False,txt='-',**default_param)
    atoms.set_calculator(calc)
    e = atoms.get_potential_energy()
    assert abs(e-6*e1)<eps
    
    
    promC = calc.get_promotion_energy(0)
    promH = calc.get_promotion_energy(1)
    AC = calc.get_atom_energy(0)
    AH = calc.get_atom_energy(1)
    BCC = calc.get_bond_energy(0,2)
    BCH = calc.get_bond_energy(0,1)
    ABC = calc.get_atom_and_bond_energy(0)
    ABH = calc.get_atom_and_bond_energy(1)
    # the energies are different here in C6H6-case than in 
    # the CH-case, because energies are divided differently
    # into bonding and atoms
    #print '\nC6H6'
    #print 'prom',promC,promH
    #print 'atom',AC, AH
    #print 'bond',BCC,BCH
    #print 'ab',ABC,ABH
    assert abs(e-sum([calc.get_atom_and_bond_energy(i) for i in range(12)] ))<eps
    
    #
    # graphene
    #
    atoms = graphene(1,1,1.42)
    calc = Hotbit(SCC=False,kpts=(10,10,1),txt='-',**default_param)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    ABC = calc.get_atom_and_bond_energy(0)
    assert abs(ABC--9.62578724411)<eps
    

  

    
if False:
    #
    #   Mayer bond order for benzene:
    # 1) model benzene with CH
    # 2) model benzene with full C6H6
    # -->>> Mayer bond-order seems not be valid with CH-modeling 
    # -->>> related to orbitals overlapping with themselves????
    #
    # THIS COULD BE INVESTIGATED MORE (9.4 2010)
    #
    #             FIXME
    #
    # using CH to model benzene
    a = Atoms('CH',[(1.395,0,0),(2.482,0,0)],container='Wedge')
    a.container.set(M=6,height=3)
    #view(a)
    c = a.extended_copy(((-2,3),1,1))
    #view(c)
    calc = Hotbit(SCC=True,kpts=(1,1,1),txt='-',**default_param)
    a.set_calculator(calc)
    e6 = a.get_potential_energy()*6
    
    MCC = calc.get_mayer_bond_order(0,1) 
    MCH = calc.get_mayer_bond_order(0,0) 
    MHH = calc.get_mayer_bond_order(1,1) 
    print(MCC, MCH, MHH)
    
    # Using the full benzene
    atoms = a.extended_copy(((-2,3),1,1))
    calc = Hotbit(SCC=True,txt='-',**default_param)
    atoms.set_calculator(calc)
    atoms.translate(-atoms.get_center_of_mass())
    e = atoms.get_potential_energy()
    assert abs(e-e6)<eps
    
    MCC = calc.get_mayer_bond_order(0,2) 
    MCH = calc.get_mayer_bond_order(0,1)
    MHH = calc.get_mayer_bond_order(1,3)
    assert abs(MCC-1.42268267342)<eps 
    assert abs(MCH-0.959404745713)<eps
    assert abs(MHH-0.00519237083333)<eps
    
    
    
    
        
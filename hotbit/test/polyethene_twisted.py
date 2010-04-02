from ase import *
from hotbit import *
from numpy import *
from box.md import check_energy_conservation
from box.fd_forces import check_forces

###

d     = 1.4
debug = False

###

atoms = Atoms('C4H4',[(0,0,0),(0.5,0,d),(0,0,2*d),(0.5,0,3*d),\
                      (-1,0,0),(1.5,0,d),(-1,0,2*d),(1.5,0,3*d)],\
#                  container='Bravais',
#              cell=[4*d,4*d,4*d], pbc=[False,False,True])
#atoms.center()
                    container='Chiral')
atoms.set_container(angle=2*pi/30,height=4*d)
#view(atoms)

calc1 = Hotbit(SCC          = False,
               kpts         = (1,1,20),
               txt          = 'polyethene.cal')
calc2 = Hotbit(SCC          = True,
               gamma_cut    = 3*d,
               kpts         = (1,1,20),
               verbose_SCC  = True,
               # For the default 1e-3 I obtain an error in force of 2 eV/A!
               mixer        = { 'name': 'anderson', 'convergence': 1e-6 },
               txt          = 'polyethene.cal')
calc3 = Hotbit(SCC             = True,
               coulomb_solver  = 'me',
               verbose_SCC     = True,
               kpts            = (1,1,20),
               mixer           = { 'name': 'anderson', 'convergence': 1e-5 },
               txt             = 'polyethene.cal')

for calc in [ calc2 ]:
    atoms.set_calculator(calc)

    # Relax (twist) the structure
    q = FIRE(atoms,trajectory='polyethene.trj',logfile=None)
    q.run(fmax=0.5)

    # Displace atoms from their equilibrium positions and check forces
    atoms.rattle(0.1)
    ffd, f0, err = check_forces(atoms)
    if debug:
        print "Finite differences forces:"
        print ffd
        print "Analytical forces:"
        print f0
        print "Error:"
        print err

    assert err < 0.01

    # Check energy conservation from a molecular dynamics run
    assert check_energy_conservation(atoms,dt=0.25*fs,steps=100,
                                     tol=0.01,plot=debug)


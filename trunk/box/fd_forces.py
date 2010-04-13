"""
Finite-differences computation of forces and the virial. Required by the
test-suite to check consistency of forces, the virial and energy.
"""

import numpy as np

import ase

def check_forces(atoms, dx=1e-6):
    """
    Compute forces and compare to forces computed numerically from a 
    finite differences approach.
    """

    f0   = atoms.get_forces().copy()
    ffd  = f0.copy()

    for a in atoms:
        r0  = a.get_position().copy()

        a.set_x(r0[0]-dx)
        ex1  = atoms.get_potential_energy()
        a.set_x(r0[0]+dx)
        ex2  = atoms.get_potential_energy()
        a.set_x(r0[0])

        a.set_y(r0[1]-dx)
        ey1  = atoms.get_potential_energy()
        a.set_y(r0[1]+dx)
        ey2  = atoms.get_potential_energy()
        a.set_y(r0[1])

        a.set_z(r0[2]-dx)
        ez1  = atoms.get_potential_energy()
        a.set_z(r0[2]+dx)
        ez2  = atoms.get_potential_energy()
        a.set_z(r0[2])

        ffd[a.index, 0]  = -(ex2-ex1)/(2*dx)
        ffd[a.index, 1]  = -(ey2-ey1)/(2*dx)
        ffd[a.index, 2]  = -(ez2-ez1)/(2*dx)

    df     = ffd-f0

    return ( ffd, f0, np.max(np.abs(df)) )


def check_field(atoms, calc, dx=1e-6):
    """
    Compute electostatic field and compare to electrostatic field computed
    numerically from a finite differences approach from the electrostatic
    potential.
    """

    E0   = calc.get_field(atoms).copy()
    Efd  = E0.copy()

    for a in atoms:
        r0  = a.get_position().copy()

        a.set_x(r0[0]-dx)
        px1  = calc.get_potential(atoms)[a.index]
        a.set_x(r0[0]+dx)
        px2  = calc.get_potential(atoms)[a.index]
        a.set_x(r0[0])

        a.set_y(r0[1]-dx)
        py1  = calc.get_potential(atoms)[a.index]
        a.set_y(r0[1]+dx)
        py2  = calc.get_potential(atoms)[a.index]
        a.set_y(r0[1])

        a.set_z(r0[2]-dx)
        pz1  = calc.get_potential(atoms)[a.index]
        a.set_z(r0[2]+dx)
        pz2  = calc.get_potential(atoms)[a.index]
        a.set_z(r0[2])

        Efd[a.index, 0]  = -(px2-px1)/(2*dx)
        Efd[a.index, 1]  = -(py2-py1)/(2*dx)
        Efd[a.index, 2]  = -(pz2-pz1)/(2*dx)

    dE = Efd-E0

    return ( Efd, E0, np.max(np.abs(dE)) )


def check_virial(atoms, de=1e-6):
    """
    Compute virial and compare to virial computed numerically from a 
    finite differences approach.
    """

    s0   = atoms.get_stress().copy()
    V0   = atoms.get_volume()
    sfd  = np.zeros([ 3, 3 ])
    c0   = atoms.get_cell().copy()

    un       = np.zeros([3,3])
    un[0,0]  = 1.0
    un[1,1]  = 1.0
    un[2,2]  = 1.0


    for i in range(3):
        for j in range(3):
            c          = c0.copy()
            eps        = un.copy()

            eps[i, j]  = un[i, j]-de
            c          = np.dot(c0, eps)
            atoms.set_cell(c, scale_atoms=True)
            e1  = atoms.get_potential_energy()

            eps[i, j]  = un[i, j]+de
            c          = np.dot(c0, eps)
            atoms.set_cell(c, scale_atoms=True)
            e2  = atoms.get_potential_energy()

            sfd[i, j]  = (e2-e1)/(2*de)

    sfd  = np.array( [ sfd[0,0], sfd[1,1], sfd[2,2], (sfd[1,2]+sfd[2,1])/2, (sfd[0,2]+sfd[2,0])/2, (sfd[0,1]+sfd[1,0])/2 ] )/V0

    return ( sfd, s0, np.max(np.abs(sfd-s0)) )


#
#    Auxiliary functions for distorting nanomaterials
#
import numpy as nu
from hotbit import Atoms


def setup_bending(atoms,angle,radius,rotation=0.0,physical=True):
    """
    Prepare a bending setup for a tube or slab.
    
    Tube should be originally periodic in z-direction with
    the correct cell, and optionally periodic in y-direction.
    
    Then atoms are set up for bending with hotbit.Wedge class
    with bending wrt. z-axis, and optionally periodic along
    z-axis.
    
    Atoms are first rotated pi/2 around -x-axis, then
    rotated an angle 'rotation' around y-axis, then
    transformed the distance 'radius' towards x-axis.
    The atoms are further adjusted for the wedge angle 'angle'.
    'physical' tells if the angle should be 2*pi/integer.
    """
    a = atoms.copy()
    a.rotate('-x',nu.pi/2)
    L = a.get_cell().diagonal()
    if abs(L.prod()-a.get_volume())>1E-10:
        raise AssertionError('Cell should be orthorhombic.')
    pbc = a.get_pbc()
    if not pbc[2] or pbc[0]:
        raise AssertionError('Should be periodic in z-direction and not periodic in x-direction')
    r = a.get_positions()    
    if any( r[:,1]<0 ):
        raise AssertionError('For bending, all atoms should be above xy-plane')
    
    # move and adjust
    a.rotate('y',rotation)
    for i in range(len(a)):
        x,y = r[i,0:2]
        R = x + radius
        phi = y*angle/L[2]
        r[i,0] = R*nu.cos(phi)
        r[i,1] = R*nu.sin(phi)
    a.set_positions(r)
    a = Atoms(a,container='Wedge')
    a.set_container(angle=angle,height=L[1],pbcz=pbc[1],physical=physical)
    return a
    

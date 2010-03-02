from ase import *
from hotbit import *
from numpy import *
from box.systems import nanotube

atoms = nanotube(9,0)
angle = atoms.container.get('angle')
traj = PickleTrajectory('twisting.traj','w',atoms)

# twist without scaling
for twist in linspace(0,pi/10,100):
    atoms.set_container(angle=angle+twist)
    traj.write()
    
# twist with scaling
atoms.set_container(angle=angle)
for twist in linspace(0,pi/10,100):
    atoms.set_container(angle=angle+twist,scale_atoms=True)
    traj.write()
    
# twist with scaling + view copies
cp = atoms.extended_copy((1,1,10))
traj = PickleTrajectory('twisting_extended.traj','w',cp)

atoms.set_container(angle=angle,scale_atoms=True)
for twist in linspace(0,pi/10,100):
    atoms.set_container(angle=angle+twist,scale_atoms=True)
    cp.set_positions( atoms.extended_copy((1,1,10)).get_positions() )
    traj.write()
    
    
    


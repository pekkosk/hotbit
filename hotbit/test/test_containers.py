from ase import *
from ase.io import Trajectory
from hotbit import *
from numpy import *
from box.systems import nanotube

# testing chiral
atoms = nanotube('C',1.42,5,0)
atoms = Atoms(atoms,container='Chiral')
atoms.set_container(angle=0.0)
traj=Trajectory('tmp.trj','w',atoms)
for angle in concatenate( (linspace(0,pi/2,50),linspace(pi/2,0,50)) ):
    atoms.set_container(angle=angle,scale_atoms=True)
    traj.write()
        
    
h = atoms.get_cell()[2,2]
for height in concatenate( (linspace(h,2*h,59),linspace(2*h,h,50)) ):
    atoms.set_container(height=height,scale_atoms=True)
    traj.write()
    
for height,angle in zip( linspace(h,2*h),linspace(0,pi/4,50) ):
    atoms.set_container(angle=angle,height=height,scale_atoms=True)
    traj.write()
    
    
# testing wedge
atoms = Atoms('C2',[(2,0,0),(2.5,1,1)],container='Wedge')
atoms.set_container(angle=pi/4,height=2.0)

traj=Trajectory('tmp2.trj','w',atoms)
for height in concatenate( (linspace(2,4,50),linspace(4,2,50)) ):
    atoms.set_container(height=height,scale_atoms=True)
    traj.write()

    
for M in concatenate( (range(100,3,-1),range(4,101)) ):
    atoms.set_container(M=M,scale_atoms=True)
    traj.write()
    
for height,M in zip( concatenate( (linspace(2,4,50),linspace(4,2,50)) ),concatenate( (range(50,3,-1),range(4,50)) )):
    atoms.set_container(height=height,M=M,scale_atoms=True)
    traj.write()
    

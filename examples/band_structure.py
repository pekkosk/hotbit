from ase import *
from hotbit import *
from ase.lattice.cubic import FaceCenteredCubic
import numpy as nu
import pylab
from box.interpolation import interpolate_path
pi = nu.pi

d = 4.08
atoms = FaceCenteredCubic('Au', latticeconstant=d,
                          directions=((0,1,1),(1,0,1),(1,1,0)),
                          align=False)
                          
calc = Hotbit(SCC=False, kpts=(8,8,8), txt='bs.cal')
atoms.set_calculator(calc)
atoms.get_potential_energy()

# reciprocal lattice vectors
V = atoms.get_volume()
a, b, c = atoms.get_cell()
a_ = 2*pi/V*nu.cross(b, c)
b_ = 2*pi/V*nu.cross(c, a)
c_ = 2*pi/V*nu.cross(a, b)

gamma = nu.array( (0,0,0) )

# a path which connects the L-points and gamma-point
kpts, distances, label_points = interpolate_path((-a_/2, a_/2, gamma, -b_/2, b_/2, gamma, -c_/2, c_/2, gamma, (-a_-b_-c_)/2, (a_+b_+c_)/2), 500)
labels = ["$-a$",'$a$','$\Gamma$','$-b$','$b$','$\Gamma$','$-c$','$c$','$\Gamma$','$-a-b-c$','$a+b+c$']


eigs = calc.get_band_energies(kpts,shift=True,rs='k')
e_min = nu.min(eigs)
e_max = nu.max(eigs)


# indices of the bands to plot [0,1,2,...]
bands = range(len(eigs[0])) # plot all bands

for j in bands:
    pylab.plot(distances, eigs[:,j])
    pylab.scatter(distances, eigs[:,j])
pylab.xticks(label_points, labels)
for k in label_points:
    pylab.plot([k,k], [e_min, e_max], linestyle=':', color='black')
pylab.show()


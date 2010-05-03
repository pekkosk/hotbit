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
#view(atoms)
a, b, c = atoms.get_cell()
V = nu.dot(a, nu.cross(b,c))

# reciprocal lattice vectors
a_ = 2*pi/V*nu.cross(b, c)
b_ = 2*pi/V*nu.cross(c, a)
c_ = 2*pi/V*nu.cross(a, b)

gamma = nu.array( (0,0,0) )

# a path which connects the L-points and gamma-point
kpts, distances, label_points = interpolate_path((-a_/2, a_/2, gamma, -b_/2, b_/2, gamma, -c_/2, c_/2, gamma, (-a_-b_-c_)/2, (a_+b_+c_)/2), 500)
labels = ["$-a$",'$a$','$\Gamma$','$-b$','$b$','$\Gamma$','$-c$','$c$','$\Gamma$','$-a-b-c$','$a+b+c$']

# kappa=False, give explicit k-points instead of kappa-points
calc = Hotbit(SCC=False, kpts=kpts, kappa=False)
atoms.set_calculator(calc)

atoms.get_potential_energy()

kpts = calc.st.k
eigs = calc.get_eigenvalues(atoms)
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


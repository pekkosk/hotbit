from ase import *
from hotbit import *
from pylab import *
from ase.io.trajectory import Trajectory


atoms = Atoms('Au2',[(0,0,0),(2.2,0,0)])
atoms.center(vacuum=10)
#traj = PickleTrajectory('dimer_curve.traj','w',atoms)
traj = Trajectory('dimer_curve.traj','w',atoms)

R = [2.2,2.4,2.54,2.8,3.0,3.2,3.4]
E = [-1.06,-2.08,-2.22,-1.99,-1.66,-1.31,-1.00]

class Calc:
    def __init__(self):
        pass
    
    def set(self,e):
        self.e = e
   
    def get_potential_energy(self,atoms):
        return self.e
    
    def get_forces(self,atoms):
        return None
    
    def get_stress(self,atoms):
        return None

calc = Calc()
atoms.set_calculator(calc)

for r,e in zip(R,E):
    atoms[1].x=atoms[0].x+r
    calc.set(e)
    print(atoms.get_potential_energy())
    traj.write()

plot(R,E)
calc = Hotbit()
atoms.set_calculator(calc)
E2=[]
R=linspace(2.3,4,50)
for r in R:
    atoms[1].x=atoms[0].x+r
    E2.append(atoms.get_potential_energy())
    
ylim(ymax=0.0)
plot(R,E2)
show()
    

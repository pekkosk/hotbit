from ase import *
from box import Atoms
from hotbit import Calculator
from hotbit import Calculator0
from ase.units import Bohr, Hartree

conv=1E-4
SCC=True
calc0=Calculator0(SCC=True,txt='test.cal',Anderson_memory=0,convergence=conv)
calc=Calculator(SCC=True,txt='test.cal',Anderson_memory=0,verbose_SCC=False,convergence=conv)

atoms=Atoms(symbols='C3H2N2H',positions=\
    [(2.323418,      1.406138,      1.072511),
    (2.250634,      2.104251,      2.171649),
    (3.297695,      3.006272,      2.785180),
    (2.879859,      4.033730,      2.950988),
    (1.315771,      2.088466,      2.823912),
    (5.310994,      3.247413,      1.364877),
    (4.437990,      3.136665,      1.974815),
    (3.628584,      2.595820,      3.775619)])
atoms.set_cell([8*Bohr,6*Bohr,6*Bohr],fix=True)
atoms.set_pbc(True)
atoms0=atoms.copy()


atoms0.set_calculator(calc0)
atoms0.get_potential_energy()
atoms.set_calculator(calc)
calc.solve_ground_state(atoms)
ddq=max(abs(calc0.get_dq()-calc.get_dq()))
de=max(abs(calc0.get_eigenvalues()-calc.get_eigenvalues()))
dn=max(abs(calc0.get_occupations()-calc.get_occupations()))
debs=abs(calc.get_band_structure_energy()-calc0.get_band_structure_energy())
decoul=abs(calc0.get_coulomb_energy()-calc.get_coulomb_energy())
de=abs(calc0.get_potential_energy(atoms)-calc.get_potential_energy(atoms))
assert ddq<1E-3
assert de<0.01
assert dn<1E-6
assert debs<1E-3
assert decoul<1E-3
assert de<1E-3




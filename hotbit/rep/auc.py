from hotbit import Calculator
from hotbit import RepulsivePotential
from box import Atoms
import numpy as nu
import ase

calc=Calculator(SCC=True,convergence=1E-9,verbose=False,\
                txt='fitting.cal',tables={'AuC':'Au_C.par'})
rep=RepulsivePotential(['C','Au'],r_cut=2.6459,r_dimer=1.85899,calc=calc)
rep.use_dimer(1.0)
rep.add_point([2.2,-8,0.8,'added'])
    
edft=[-1533.0049,-1538.3640,-1539.8517,-1539.9728,-1539.8882,-1539.4188,-1538.5373]

d=1.42    
h=d*nu.cos(nu.pi/6)
rl=[1.3809,1.5809,1.7809,1.8809,1.9809,2.1809,2.4809]

traj=[]
for r in rl:
    atoms=Atoms('C6Au',[(d/2,h,0),(-d/2,h,0),(d/2,-h,0),(-d/2,-h,0),\
                        (d,0,0),(-d,0,0),(0,0,nu.sqrt(r**2-d**2))],cell=(10,10,10))
    traj.append(atoms)
    ase.view(atoms)
    print r,atoms.get_positions()
    raw_input('dsds')
    
assert False
         
#rep.use_energy_curve('traj.xyz',edft,charge=0,bonds=[(0,1)])
#rep.use_equilibrium_structure('cluster.xyz')
rep.fit()
rep.write_to_par('koe.par')
rep.plot()
#rep.write_to_file()


    
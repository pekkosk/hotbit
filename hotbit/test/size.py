from hotbit import Calculator
from hotbit import Calculator0
from ase import *
from hotbit.test.misc import default_param
from box.mix import Timer
          
from ase.lattice.cubic import FaceCenteredCubic

t=Timer()

#atoms = FaceCenteredCubic(directions=[[1,-1,0], [1,1,-2], [1,1,1]],
                          #size=(1,1,1), symbol='Au', pbc=(1,1,1))
atoms=read('zz8H.xyz')
    
atoms.set_pbc((True,False,False))
atoms.set_cell((4*2.45951,1000,3000))
atoms.center()
#view(atoms)                      
#assert False
                          
default_param['Anderson_memory']=0
default_param['convergence']=1E-2
default_param['mixing_constant']=0.1
default_param['width']=0.1
#calc0=Calculator0(verbose=True,SCC=True,gamma_cut=3,txt='sizef.cal',**default_param)
#atoms.set_calculator(calc0)
##atoms.get_forces()
#t()
#calc0.finalize()
                      
calc=Calculator(verbose=True,SCC=True,gamma_cut=3,verbose_SCC=True,txt='size.cal',**default_param)
atoms.set_calculator(calc)
atoms.get_forces()

calc.__del__() 
#print calc.st.solver.get_iteration_info()
t()



#systems=['H2COH','AuC']          

#default_param['convergence']=1E-5
#default_param['Anderson_memory']=3
#default_param['width']=0.01
  
  
#print '    ... forces for %s, SCC=' %system, SCC         

#atoms=molecule(system)
#atoms.center(vacuum=5)
#atoms[0].z+=0.2
#atoms=Atoms(atoms)
#atoms.set_calculator(calc)

#rec, de, du=energy_conservation(atoms,dt=0.2*fs,steps=50)
#calc.__del__()
#print de/du
#assert de/du<0.01

    

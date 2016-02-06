from ase import *
from pylab import *
from hotbit import *
from box import mix
from hotbit.analysis import LinearResponse

# equilibrium length
Na2=Atoms('Na2',[(0,0,0),(3.0,0,0)],cell=(6,4,4),pbc=False)
Na2.center()
calc=Hotbit(parameters=testpar,SCC=True,txt='optical.cal')
Na2.set_calculator(calc)
calc.solve_ground_state(Na2)
lr=LinearResponse(calc,energy_cut=10,txt='lr.out')
lr.run()
lr.info()


#get excitation energies and oscillator strengths
omega, F=lr.get_linear_response() 
e,f=mix.broaden(omega,F,width=0.1,function='lorentzian')
plot(e,f)
savefig('lr.pdf')


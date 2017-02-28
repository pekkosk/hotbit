from ase import *
from hotbit import Hotbit
from ase.lattice.cubic import FaceCenteredCubic
from hotbit.test.misc import default_param
import pylab as pl
from box import Atoms



for SCC in [True,False]:
    timings=[]
    norbs=[]
    sccs=['non-SCC','SCC'][SCC]
    for nx in [1,2,3,4,5]:
        for ny in [1,2,3,4,5]:
            print('nx',nx,ny)
            atoms = FaceCenteredCubic(directions=[[1,-1,0], [1,1,-2], [1,1,1]],\
                                    size=(nx,ny,1), symbol='Au', pbc=(0,0,0))
                                    
            calc=Hotbit(verbose=True,SCC=SCC,gamma_cut=3,txt='size2.cal',**default_param)                
            atoms.set_calculator(calc)
            try:
                atoms.get_forces()
            except:
                continue                
            calc.timer.summary()
            timings.append( calc.timer.get_timings() )
            norbs.append( calc.el.get_nr_orbitals() )
            calc.__del__()
        
    order=argsort(norbs)
    norbs=sort(norbs)
    timings2=[]
    for i in order:
        timings2.append(timings[i])
    
    times={}
    maxtime=0.0
    for key in timings2[0]:    
        times[key]=array([ts[key] for ts in timings2])
        maxtime=max(maxtime,times[key].max())
        
    s=0.1        
    for key in timings[0]:
        s+=0.1
        if times[key].max()<0.05*maxtime: continue        
        pl.plot(norbs,times[key],label=key,lw=s)
        
    atoms=Atoms(atoms)    
    pl.title('Timings up to %s with %s' %(atoms.get_name(),sccs) )
    pl.xlabel('Number of orbitals')
    pl.ylabel('Time (s)')
    pl.legend(loc='upper left')    
    pl.savefig('size_timing_%s.png' %sccs)
    #pl.plot()   
    #pl.show()
    pl.clf()
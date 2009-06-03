#def quench(atoms,name=None,fmax=0.05,calc={}):
    #"""
    #Quench atoms and write quenching into traj-files and quenched .xyz and .gpw.
    
    #Return e [atoms] [calc] depending whether atoms and calc are new objects.
    #"""
    #import ase
    #from hotbit import Calculator
    #from box import mix
    
    #if type(atoms)==type(''):
        #mode='file'
    #else:
        #mode='atoms'
        
    #if mode=='file':
        #name,ext=mix.base_and_extension(atoms)
        #atoms=ase.read(atoms)
        
    #if type(calc)==type({}):
        #cmode='new'
        #kwargs={}
        #kwargs.update(calc)
        #hb=Calculator(None,**kwargs)
    #else:
        #cmode='old'
        #hb=calc
    
    #atoms.set_calculator(hb)
    #traj=ase.PickleTrajectory(name+'_quenching.trj','w',atoms)
    
    #opt=ase.QuasiNewton(atoms)
    ##opt=ase.FIRE(atoms)
    #opt.attach(traj.write)   
    #opt.run(fmax=fmax)
    #e=atoms.get_potential_energy()
    #ase.write(name+'_quenched.xyz',atoms)
    
    #if mode=='file' and cmode=='new':
        #return e,atoms,gpw
    #elif mode=='file' and cmode=='old':
        #return e,atoms
    #elif mode=='atoms' and cmode=='new':
        #return e,gpw
    #elif mode=='atoms' and cmode=='old':
        #return e  

def dimer_curve(symbols,a,b,calc={},view=False):
    """
    Calculate the dimer curve [a,b] for given atoms, and, 
    optionally, view it.
    """
    from ase import Atoms
    atoms=Atoms(symbols=symbols)
    atoms.set_calculator(calc)
    e=[]
    rlist=linspace(a,b,50)
    for r in rlist:
        atoms.set_positions([(0,0,0),(r,0,0)])
        e.append( atoms.get_potential_energy() )
        
    if view:
        import pylab as p
        p.plot(rlist,e)
        p.show()
    return rlist,e

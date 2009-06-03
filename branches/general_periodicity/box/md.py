import ase
import numpy as nu
from box import mix
from box.data import data
from ase import VelocityVerlet
from ase import PickleTrajectory
from ase.units import fs, Hartree

class TrajectoryRecording:
    def __init__(self,atoms,verbose=False):
        self.atoms=atoms
        self.properties={'epot':self.atoms.get_potential_energy,\
                         'etot':self.atoms.get_total_energy   }
        self.data={}
        for key in self.properties:
            self.data[key]=[]
        self.i=0
        self.verbose=verbose
        
    def __call__(self):
        for key in self.properties:
            self.data[key].append( self.properties[key]() )       
        self.i+=1
        if self.verbose:
            print 'record',self.i,self.data['epot'][-1]
        
    def get(self,key):
        return self.data[key]

    def plot_energies(self):
        pl.plot( self.get('epot'),label='potential')
        pl.plot( self.get('etot'),label='total')
        pl.show()
        pl.close()
        
        
def dimer_curve(atoms,r1=None,r2=None,N=100):
    """ Calculate dimer curve """
    rc=[]
    for s in atoms.get_chemical_symbols():
        rc.append(data[s]['R_cov'])
    if r1==None:
        r1=(rc[0]+rc[1])*0.5
    if r2==None:
        r2=(rc[0]+rc[1])*2.5              
    e=[]
    rl=nu.linspace(r1,r2,N)
    for r in rl:
        atoms.set_positions([(0,0,0),(r,0,0)])
        e.append( atoms.get_potential_energy() )
    return rl,nu.array(e)                    
                    
    

def canonical(atoms,dt=1.0,steps=1000,output=10,name=None,verbose=False):
    """ Perform little canonical simulation. 
    
    parameters:
    -----------
    atoms: atoms with calculator attached
    dt: time step in fs
    steps: how many md steps
    output: output frequency
    name: TrajectoryRecording name
    verbose: increase verbosity    
    """
    if name==None:
        try:
            name=atoms.get_name()
        except:
            name='microcanonical'
    name+='.trj'        
    traj=PickleTrajectory(name,'w',atoms)
    rec=TrajectoryRecording(atoms,verbose)
    md=VelocityVerlet(atoms,dt*fs)
    md.attach(rec,interval=output)  
    md.attach(traj.write,interval=output)  
    md.run(steps) 
    return rec
    

def microcanonical(atoms,dt=1.0,steps=100,output=1,name=None,verbose=False):
    """ Perform little microcanonical simulation. 
    
    parameters:
    -----------
    atoms:
    dt: time step in fs
    steps: how many md steps
    output: output frequency
    name: TrajectoryRecording name
    verbose: increase verbosity
    
    Return TrajectoryRecording object for further analysis.
    """
    if name==None:
        try:
            name=atoms.get_name()
        except:
            name='microcanonical'
    name+='.trj'        
    traj=PickleTrajectory(name,'w',atoms)
    rec=TrajectoryRecording(atoms,verbose)
    md=VelocityVerlet(atoms,dt*fs)
    md.attach(rec,interval=output)
    md.attach(traj.write,interval=output)    
    md.run(steps)    
    return rec
    
    
def energy_conservation(atoms,dt=1.0,steps=100,name=None,verbose=False):
    """ Make microcanonical simulation, check energy conservation. 
    
    parameters:
    -----------
    atoms: atoms object with calculator attached (and initial geometry)
    dt: time step to perform md (fs)
    steps: how many md steps are taken
    name: name for the TrajectoryRecording object
    verbose: increase verbosity
    
    Return TrajectoryRecording instance and deviations in etot and epot.    
    """
    rec=microcanonical(atoms,dt,steps=steps,output=1,name=name,verbose=verbose)
    etot, epot=rec.get('etot'), rec.get('epot')
    de=nu.sqrt( nu.var(etot) )
    du=nu.sqrt( nu.var(epot) )
    return rec, de, du
    
    
def energy_conservation_plot(atoms,dts,total_time=100):
    """ Make set of NVE simulations with different time steps and plot errors. 
    
    parameters:
    -----------
    atoms: atoms object with calculator attached (and initial geometry)
           perform all simulations from same initial point
    dts: list of time steps to try (fs)
    total_time: total simulation time (in fs) for all simulations
                (smaller time stepping makes more calculations)
    """
    import pylab as pl
    R0=atoms.get_positions()
    dE=[]
    dU=[]
    for dt in dts:
        steps=int(total_time/dt)
        print 'time step %.2f fs, steps %i' %(dt,steps)
        atoms.set_positions(R0)
        rec, de, du=energy_conservation(atoms,dt,steps=steps)
        dE.append(de)
        dU.append(du)
    atoms.set_positions(R0)        
    frac=nu.array(dE)/nu.array(dU)        
    rec.plot_energies()
    
    pl.plot(dts,frac*100)
    pl.title('Total energy deviation/potential energy deviation')
    pl.xlabel('time step (fs)')  
    pl.ylabel('dev(e)/dev(u) (%)')
    pl.plot()          
    pl.show()
    
 
class TrajectoryWriter:
    """
    Trajectory writer for transition paths. 
    Write each path into its own file.
    """
    def __init__(self,images,name=None):
        """ images is the list of Atoms objects; you can supply also the name. """
        
        if name==None:
            self.name='path'
        else:
            self.name=name
        self.images=images
        self.i=0
        
    def __call__(self):
        """ Writes trajectory file for current atoms list. """
        from ase import PickleTrajectory
        traj = PickleTrajectory('%s_it%i.trj' %(self.name,self.i), 'w')
        for image in self.images:
            traj.write(image)
        self.i+=1 
        
            
def quench(atoms,name=None,fmax=0.05,method='QN'):
    """ Quench given Atoms object with given attached calculator. """
    if name==None:
        try:
            name=atoms.get_name()
        except:
            raise ValueError('name not specified')
    traj=ase.PickleTrajectory(name+'_quenching.trj','w',atoms)
    
    if method=='QN':
        qn=ase.QuasiNewton(atoms)
    elif method=='FIRE':
        qn=ase.FIRE(atoms)        
    qn.attach(traj.write)   
    qn.run(fmax=fmax)
    e=atoms.get_potential_energy()
    ase.write(name+'_quenched.xyz',atoms)
    return e
    
    
def transition_barrier(calc,quench,guess,cell=None,pbc=None,constraints=None,\
                       M=None,name=None,fmax=0.05,steps=1000000):
    import ase
    method=ase.NEB
    
    if type(guess)!=type([]) and type(guess)!=type(''):
        # all Atoms properties should be set for Atoms-type guess
        images=guess
        path=method(images)
    else:
        if type(guess)==type(''):
            assert guess.split('.')[-1]=='trj' and M==None 
            # for some reason copying has to be done...
            images=[]
            for image in ase.PickleTrajectory(guess):
                images.append(image.copy())
            path=method(images)
        else:
            assert type(guess)==type([]) and M!=None
            images=[]
            first=ase.read(guess[0])
            images.append(first)
            for i in range(M-2):
                images.append(first.copy())
            last=ase.read(guess[1])
            images.append(last)
            path=method(images)
            path.interpolate()  
              
        # now coordinates are set; set other properties
        assert cell!=None and pbc!=None 
        for image in images:
            image.set_pbc(pbc)
            image.set_cell(cell,fix=True)
            if constraints!=None:
                image.set_constraint(constraints)
                
    # attach calculators
    if type(calc)==type([]): 
        for image,clc in zip(images,calc):
            image.set_calculator(clc)   
    else:    
        for image in images:
            image.set_calculator(calc.copy())   
            
    for image in images:
        image.center()            
            
    if quench:
        e1=quench_atoms(images[0],'initial')
        print 'First image quenched. Energy=',e1
        eM=quench_atoms(images[-1],'final')
        print 'Last image quenched. Energy=',eM
            
    # solve the transition path
    writer=TrajectoryWriter(images)        
    minimizer=ase.QuasiNewton(path)
    minimizer.attach(writer)
    minimizer.run(fmax=fmax,steps=steps)
    
    # output of the final path
    traj = ase.PickleTrajectory('%s_converged.trj' %name, 'w')
    path.write(traj)
    return images
       
def merge_atoms(atoms1,atoms2,box_from=None,R=0.3):
    """ 
    Merge two sets of atoms without atom duplication.
    (atoms1 and atoms2 may share same atoms within radius R)
    """
    from box import Atoms
       
    atoms=atoms1.copy()
    print type(atoms)

    # set box from the first atoms or box_from
    if box_from!=None:
        bx=box_from
    else:
        bx=atoms1 
    atoms.set_cell(bx.get_cell(),fix=True)
    atoms.set_pbc(bx.get_pbc())
    
    # add atoms from atoms2 if not close
    for a2,r2 in zip(atoms2,atoms2.get_positions()):
        close=False
        for r1 in atoms1.get_positions():
            if atoms.distance(r1,r2)<R:
                close=True
        if not close:
            atoms=atoms+a2    
    return atoms
    
    


if __name__=='__main__':
    from box import Atoms
    atoms1=Atoms(symbols='C4H',positions=[(0,0,0),(0,0,1),(0,0,2),(0,0,3),(0,0,4)],cell=(7,8,9),pbc=True)
    atoms2=Atoms(symbols='C4H',positions=[(0,0,0.29),(0,0,1),(0,0,7),(0,0,8),(0,0,9)],cell=(170,18,19),pbc=True)
    atoms=merge_atoms(atoms1,atoms2,box_from=atoms2)
    print atoms.get_positions()
    print atoms.get_cell()
    print atoms.get_pbc()
    
            
            
            
            
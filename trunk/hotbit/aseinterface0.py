"""
    ASE-calculator interface for HOTBIT.
    
    P. Koskinen 21.1 2008
    
"""

import os
import sys
import weakref

import numpy as npy
import ase
from box.timing import Timer
from ase.units import Bohr, Hartree
from elements import Elements
from interactions import Interactions
from hotbit.output import Output
import hotbit.auxil as aux
import box.mix as mix
vec=npy.array
err=mix.error_exit
Ha=27.2113956
a0=0.529177
f0=Ha/a0
import time


         
    
class Calculator0(Output):
    """
    ASE-calculator frontend for HOTBIT calculations.
    """
    def __init__(self,command=None,elements=None,tables=None,verbose=True,**kwargs):
        """ 
        Initialize calculator. 
        
        key               default
        ==========================================================================
        propagation       BO        propagation mode
        charge            0.0       total electric charge for system
        width             0.0005    width of Fermi occupation 
        SCC               True      Self-Consistent Charge?
        convergence       1E-3      convergence criterion 
                                    (Mulliken charge variation/atom)
        mixing_constant   0.2       density mixing for Anderson mixer
        electric_field    '0 0 0'   Electric field vector (as string)
        Anderson_memory   3         Memory for Anderson mixing
        maxiter           50        Maximum number of self-consistent 
                                    iterations (if SCC)
        gamma_cut         1E10      Range for Coulomb interaction
        """
        
        import os
        from copy import copy
        defaults={'propagation':'BO',\
                  'charge'            :0.0,\
                  'width'             :0.02,\
                  'SCC'               :True,\
                  'convergence'       :1E-3,\
                  'mixing_constant'   :0.2,\
                  'electric_field'    :'0 0 0',\
                  'Anderson_memory'   :3,\
                  'maxiter'           :50,\
                  'gamma_cut'         :1E10,\
                  'txt'               :None} 
            
        self.args=copy(defaults)    
        self.args.update(kwargs)
        self.args['width']=self.args['width']/Hartree
        self.command=command
        self.set_up=False
        self.element_files=elements
        self.table_files=tables    
        self.verbose=verbose     
        self.version='3.0 fortran-based'   
        
    def set_text(self,txt):
        """ Set up the output file. """
        if txt is None:
            self.txt=sys.stdout
        else:
            self.txt=open(out,'a')
        
    def set(self,**kwargs):
        """ (Re-)Set calculator parameters, _before_ calculation. """
        for arg in kwargs:
            value=kwargs[arg]
            if arg in self.args:
                self.args.update({arg:value})
            elif arg is 'elements':
                if self.element_files is None:
                    self.element_files=value
                else:
                    self.element_files.update(value)
            elif arg is 'tables':
                if self.table_files is None:
                    self.table_files=value
                else:                    
                    self.table_files.update(value)
            elif arg is 'verbose':
                self.verbose=value
            elif arg is 'command':
                self.command=value
        
    def finalize(self):
        """ Close the communication pipes for hotbit. """
        if self.set_up:
            self.read_until('hotbit waiting for positions')
            self.inp.write('final comm\n') 
            self.inp.flush()
            print>>self.txt, 'Hotbit pipes explicitly closed.'
            self.read_until('End communication & stop.')
            
            
    def read_until(self,string):
        while 1:
            line=self.txt.readline() #'' only at EOF 
            if line.find(string)>=0 or line=='': break
            if self.verbose: print>>self.txt,  '>',line,
            #print>>self.txt,  '>>>',string
    
    def get_arguments(self):
        return self.args
    
    def get(self,arg=None):
        if arg==None:
            return self.args
        else:
            return self.args[arg]
            
    def __del__(self):
        self.finalize()
      
    def _set_up_calc(self,atoms,clean=False):
        """ Make atoms.dat file for setting up hotbit. """
        self.set_text(self.args['txt'])
        self.timer=Timer('hotbit fortran',txt=self.txt)
        from os import popen2
        if self.command==None:
            self.command=os.environ.get('HOTBIT_EXE')
        self.el=Elements(atoms,self.element_files,charge=self.args['charge'])
        self.ia=Interactions(self,self.timer,self.el,self.table_files)
        self.write_md_dat(self.args,elements=self.el.get_files(),\
                          tables=self.ia.get_files())                          
        write_atoms_dat(atoms)
        self.greetings()  
                
        print>>self.txt,  'Open hotbit pipes for communication...'
        self.inp,self.txt=popen2(self.command,'t')  
        self.set_up=True
        self.atoms=atoms
        if clean:
            from time import sleep
            sleep(1.0)
            os.remove('md.dat')
            os.remove('atoms.dat')
            
    def greetings(self):
        """ Simple greetings text """
        import time
        import os
        print>>self.txt,  ' _           _    _     _ _'
        print>>self.txt,  '| |__   ___ | |_ | |__ |_| |_'
        print>>self.txt,  '|  _ \ / _ \|  _||  _ \| |  _|'
        print>>self.txt,  '| | | | ( ) | |_ | ( ) | | |_'
        print>>self.txt,  '|_| |_|\___/ \__|\____/|_|\__|',self.version
        print>>self.txt,  'Date:',time.asctime()
        dat=os.uname()
        print>>self.txt,  'Nodename:',dat[1]
        print>>self.txt,  'Arch:',dat[4]
        print>>self.txt,  'Dir:',os.path.abspath(os.path.curdir)
        print>>self.txt, self.el.greetings()
        print>>self.txt, self.ia.greetings()
        
    def copy(self):
        """ Return simple copy of this calculator. """
        from copy import deepcopy
        return deepcopy(self)
        
    def write_md_dat(self,inpt,elements,tables):
        """ Write input file to start hotbit communication. """
        f=open('md.dat','w')
        keys=inpt.keys()
        keys.sort(lambda a,b: cmp(a.lower(),b.lower()))
        for key in keys:
            print>>f, key,'=',inpt[key]
        for element in elements:
            print>>f, element,'=',elements[element]
        for table in tables:
            print>>f, table,'=',tables[table]
        f.close()    
        
    def get_potential_energy(self,atoms):
        """ Return the potential energy of present system. """
        if not self.set_up:
            self._set_up_calc(atoms)
        self.hotbit_feed(atoms)        
        self.hotbit_fetch(atoms)
        return self.energy*Ha
              
    def get_forces(self,atoms):
        """ Return the forces of present system. """
        if not self.set_up:
            self._set_up_calc(atoms)
        self.hotbit_feed(atoms)
        self.hotbit_fetch(atoms)
        return vec(self.forces)*f0
            
    def get_stress(self,atoms):
        return None       
         
    def get_charge(self):
        return None
    
    def get_dq(self):
        return self.dq
    
    def get_eigenvalues(self):
        return self.eigenvalues*Ha

    def get_occupations(self):
        return self.occupations
    
    def get_band_structure_energy(self):
        return self.ebs*Hartree
    
    def get_coulomb_energy(self):
        return self.ecoul*Hartree
   
    def get_optical(self):
        op=self.optical.copy()
        op[:,0]=op[:,0]*Hartree
        return op
            
    def hotbit_feed(self,atoms):
        self.read_until('hotbit waiting for positions')
        print>>self.inp, 'positions to hotbit'
        mix.write(self.inp,atoms.get_positions()/a0)
        mix.write(self.inp,atoms.get_cell()/a0)
        self.inp.flush()
    
    def hotbit_fetch(self,atoms):
        self.read_until('hotbit > standard feed')
        while 1:
            line=self.txt.readline().split('=')
            key=line[0].strip()
            if len(line)>1:
                value=line[1]
            if key=='epot':
                self.energy=float(value)
            elif key=='erep':
                self.erep=float(value)
            elif key=='ebs':
                self.ebs=float(value)
            elif key=='ecoul':
                self.ecoul=float(value)
            elif key=='forces':
                self.forces=[]
                for atom in atoms:
                    f1=vec([float(fi) for fi in self.txt.readline().split()])
                    self.forces.append(vec(f1))
            elif key=='dq':
                self.dq=[]
                for atom in atoms:
                    self.dq.append(float(self.txt.readline()))
                self.dq=npy.array(self.dq)
            elif key=='eigenvalues':
                self.eigenvalues=[]
                for i in range(self.el.get_nr_orbitals()):
                    self.eigenvalues.append(float(self.txt.readline()))
                self.eigenvalues=npy.array(self.eigenvalues)                
            elif key=='occupations':
                self.occupations=[]
                for i in range(self.el.get_nr_orbitals()):
                    self.occupations.append(float(self.txt.readline()))
                self.occupations=npy.array(self.occupations)
            #elif key=='optical':
                #self.optical=[]
                #while 1:
                    #line=self.txt.readline()
                    #if line.find('end optical')>=0: break
                    #x,f=line.split()
                    #self.optical.append([float(x),float(f)])
                #self.optical=npy.array(self.optical)
            elif key=='end comm':
                break
            else:
                raise RuntimeError('not valid line %s from hotbit' %line)
        self.inp.write('no additional data\n')
        self.inp.flush()
       
    def get_repulsive_energy(self,include=None):
        """ Calculate the repulsive energy with included element pairs. """
        erep=0.0     
        #print self.el.get_atom_pairs(3,include)
        for i,ei,j,ej,d in self.el.get_atom_pairs(3,include):
            if i==j: 
                continue
            erep+=self.ia.vrep[ei+ej](d)
        return erep*Hartree      
       
    def get_repulsive_forces(self,include=None):
        """ Calculate the forces due to repulsive potentials for element pairs. """
        f=npy.zeros((self.el.N,3))
        for i,ei,j,ej,d,r,rhat in self.el.get_atom_pairs(5,include):
            if i==j: 
                continue
            frep=self.ia.vrep[ei+ej](d,der=1)*rhat
            f[i,:]=f[i,:]+frep
            f[j,:]=f[j,:]-frep
        return f*(Hartree/Bohr)
        
         
        
def write_atoms_dat(atoms):
    """ Write the atoms into atoms.dat-file. """
    f=open('atoms.dat','w')
    print>>f, 'nat=',len(atoms)
    print>>f, '\natoms='
    symb=atoms.get_chemical_symbols()
    for s in symb:
        print>>f, s,0,0,0
    
    print>>f, '\n\ncell='
    cell=atoms.get_cell()
    pbc=atoms.get_pbc()
    #for i in range(3):                               
        #if cell[i,i]<1.1 or not pbc[i]: cell[i,i]=100
    #atoms.set_cell(cell,fix=True)  # do not touch Atoms!!!
    for i in range(3):
        if not pbc[i]: cell[i,i]=1E6
    cell=cell/a0
    mix.write(f,cell)        
    f.close()
    

        
        
if __name__=='__main__':
    from ase import *
    from box import Atoms
    #h2o=Atoms(symbols='C3H2N2H',positions=[(0,0,0),(2,2,2),(3,3,3),(3,4,5),(1,2,3),(1,2,4),(4,3,2),(4,3,6)],cell=(10,10,10))
    calc=Calculator(SCC=True,txt='test.cal',Anderson_memory=0)
    #h2o.set_calculator(calc)
    #h2o.get_potential_energy()
    #calc.states.solve()
    
    #t=mix.Timer()
    #large=Atoms(symbols='100C',positions=[(i*1.42,0,0) for i in range(100)], cell=(1000,20,20))
    
    
    #large.set_calculator(calc)
    ##t()
    ##for i in range(10):
    #large.get_potential_energy()
    ##t()
    ##for i in range(10):
    
    #calc.states.solve()
    #t()
    atoms=Atoms(symbols='C3H2N2H',positions=\
       [(2.323418,      1.406138,      1.072511),
       (2.250634,      2.104251,      2.171649),
       (3.297695,      3.006272,      2.785180),
       (2.879859,      4.033730,      2.950988),
       (1.315771,      2.088466,      2.823912),
       (5.310994,      3.247413,      1.364877),
       (4.437990,      3.136665,      1.974815),
       (3.628584,      2.595820,      3.775619)])
    atoms.set_cell([8*a0,6*a0,6*a0],fix=True)
    atoms.set_pbc(True)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    calc.states.solve()
    
    
 
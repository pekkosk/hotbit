# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

"""
    ASE-calculator interface for HOTBIT.
    (Handy Open-Source Tight-Binding Tool)
    
"""
import weakref

import numpy as nu
import ase
from ase.units import Bohr, Hartree
from box.timing import Timer
from elements import Elements
from interactions import Interactions
from electrostatics import Electrostatics
from repulsion import Repulsion
from states import States
from hotbit.output import Output
import hotbit.auxil as aux
import box.mix as mix

vec=nu.array
err=mix.error_exit
Ha=27.2113956
a0=0.529177
f0=Ha/a0


         
    
class Calculator(Output):
    """
    ASE-calculator frontend for HOTBIT calculations.
    """
    def __init__(self,elements=None,tables=None,verbose=True,**kwargs):
        """ 
        Initialize calculator. 
        
        Parameters:
        -----------
        charge            total electric charge for system (-1 means an additional electron)
        width             width of Fermi occupation (in eV)
        SCC               Self-Consistent Charge calculation
        convergence       convergence criterion (Mulliken charge variation/atom)
        mixing_constant   density mixing for Anderson mixer
        Anderson_memory   Memory for Anderson mixing
        maxiter           Maximum number of self-consistent iterations (for SCC)
        gamma_cut         Range for Coulomb interaction
        elements          Dictionary for elements used, e.g. {'H':'H_custom.elm','C':'/../C.elm'}
                          Items can also be elements directly: {'H':H} (with H being type Element)
                          * If elements==None, use default element info (from HOTBIT_PARAMETERS).
                          * If elements['others']=='default', use default parameters for all other
                            elements than the ones specified. E.g. {'H':'H.elm','others':'default'}
                            (otherwise all elements present have to be specified excplicitly).
        tables:           Dictionary for interactions, e.g. {'CH':'C_H.par','CC':'C_C.par'}
                          * If elements==None, use default interactions.
                          * If elements['others']='default', use default parameters for all other
                            interactions than the ones specified. 
                            E.g. {'CH':'C_H.par','others':'default'}
        """       
        from copy import copy
        defaults={'charge'            :0.0,\
                  'width'             :0.02,\
                  'SCC'               :True,\
                  'convergence'       :1E-3,\
                  'mixing_constant'   :0.2,\
                  'Anderson_memory'   :3,\
                  'maxiter'           :50,\
                  'gamma_cut'         :None,\
                  'txt'               :None,\
                  'verbose_SCC'       :False} 
                  
        for key in kwargs:
            if key not in defaults:
                raise AssertionError('Not valid keyword argument %s.' %key)           
            
        self.args=copy(defaults)    
        self.args.update(kwargs)
        self.args['width']=self.args['width']/Hartree
        if self.args['gamma_cut']!=None:
            self.args['gamma_cut']=self.args['gamma_cut']/Bohr
        self.init=False
        self.element_files=elements
        self.table_files=tables    
        self.verbose=verbose     
        self.set_enabled=True
        self.version='3.0 alpha' 
        self.notes=[]
        self.set_text(self.args['txt'])
        self.timer=Timer('Hotbit',txt=self.get_output())
       
        
    def __del__(self):
        """ Delete calculator -> timing summary. """
        if self.get('SCC'):
            try:
                print>>self.txt, self.st.solver.get_iteration_info()
                self.txt.flush()
            except:
                pass
        if len(self.notes)>0:            
            print>>self.txt, 'Notes and warnings:'            
            for note in self.notes:
                print>>self.txt, note         
        self.timer.summary()
        self.close_output()
        
    def add_note(self,note):
        """ Add warning (etc) note to be printed in log file end. """
        self.notes.append(note)        
          
    def greetings(self):
        """ Simple greetings text """
        from time import asctime
        from os import uname
        from os.path import abspath, curdir
        print>>self.txt,  '\n\n\n\n\n'
        print>>self.txt,  ' _           _    _     _ _'
        print>>self.txt,  '| |__   ___ | |_ | |__ |_| |_'
        print>>self.txt,  '|  _ \ / _ \|  _||  _ \| |  _|'
        print>>self.txt,  '| | | | ( ) | |_ | ( ) | | |_'
        print>>self.txt,  '|_| |_|\___/ \__|\____/|_|\__|',self.version
        print>>self.txt,  'Date:',asctime()
        dat=uname()
        print>>self.txt,  'Nodename:',dat[1]
        print>>self.txt,  'Arch:',dat[4]
        print>>self.txt,  'Dir:',abspath(curdir)
        print>>self.txt,  'System:',self.el.get_name()
        print>>self.txt,  '       Charge=%4.1f' %self.args['charge']
        print>>self.txt,  '       Box: (Ang)', nu.array(self.el.get_box_lengths())*Bohr
        print>>self.txt,  '       PBC:',self.pbc
        print>>self.txt, self.el.greetings()
        print>>self.txt, self.ia.greetings()  
        print>>self.txt, self.rep.greetings()   
             
    def out(self,text):
        print>>self.txt, text
        
    
    def set_text(self,txt):
        """ Set up the output file. """
        if txt is None:
            from sys import stdout
            self.txt=stdout
        else:
            self.txt=open(txt,'a')
        
    def set(self,**kwargs):
        """ (Re-)Set calculator parameters, _before_ calculation. """
        if not self.set_enabled:
            raise RuntimeError('Calculator initialized -> set method disabled.')
        for arg in kwargs:
            value=kwargs[arg]
            if arg in self.args:
                if arg is 'width':
                    self.args.update({arg:value/Hartree})
                elif arg is 'gamma_cut':
                    self.args.update({arg:value/Bohr})
                else:
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
       
    def get(self,arg=None):
        if arg==None:
            return self.args
        else:
            return self.args[arg]
        
        
    def solve_ground_state(self,atoms):
        """ If atoms moved, solve electronic structure. """
        if not self.init:
            self._initialize(atoms)
        if not self.el.is_solved(atoms):
            self.st.solve()
        else:
            pass
                  
                  
    def _initialize(self,atoms):
        """ Initialization of hotbit. """  
        self.timer.start('initialization')
        self.el=Elements(self,atoms,self.timer,self.element_files,charge=self.args['charge'])
        self.ia=Interactions(self,self.timer,self.el,self.table_files)
        self.es=Electrostatics(self,self.timer)
        self.st=States(self,self.timer,self.el,self.ia)
        self.rep=Repulsion(self,self.timer,self.el,self.ia)
        self.init=True
        self.set_enabled=False
        self.pbc=atoms.get_pbc()
        self.greetings()
        if self.get('SCC') and nu.any(self.pbc) and self.get('gamma_cut')==None:
            raise NotImplementedError('SCC not implemented for periodic systems yet (see parameter gamma_cut).')
        if nu.any(self.pbc) and abs(self.get('charge'))>0.0:
            raise AssertionError('Charged system cannot be periodic.')
        self.flush()
        self.timer.stop('initialization')
        
        
    def get_potential_energy(self,atoms):
        """ Return the potential energy of present system. """
        self.solve_ground_state(atoms)
        self.timer.start('energies')
        ebs=self.get_band_structure_energy(atoms)
        ecoul=self.get_coulomb_energy(atoms)
        erep=self.rep.get_repulsive_energy()
        self.timer.stop('energies')
        return ebs+ecoul+erep
              
              
    def get_forces(self,atoms):
        """ Return the forces of present system. """
        self.solve_ground_state(atoms)
        self.timer.start('forces')
        fbs=self.st.band_structure_forces()
        frep=self.rep.get_repulsive_forces()
        fcoul=self.es.gamma_forces() #zero for non-SCC
        self.timer.stop('forces')
        return (fbs+frep+fcoul)*(Hartree/Bohr)
            
            
    def get_stress(self,atoms):
        self.solve_ground_state(atoms)
        return None            
         
         
    def get_charge(self):
        return self.get('charge')
    
    
    def get_dq(self,atoms):
        self.solve_ground_state(atoms)
        return self.st.get_dq()
    
    
    def get_eigenvalues(self,atoms):
        self.solve_ground_state(atoms)
        return self.st.get_eigenvalues()*Hartree
    
    
    def get_occupations(self):
        #self.solve_ground_state(atoms)
        return self.st.get_occupations()
        
        
    def get_band_structure_energy(self,atoms):
        self.solve_ground_state(atoms)
        return self.st.band_structure_energy()*Hartree
            
            
    def get_coulomb_energy(self,atoms):
        self.solve_ground_state(atoms)
        return self.es.coulomb_energy()*Hartree
    
    
    # some not implemented ASE-assumed methods
    def get_fermi_level(self):
        raise NotImplementedError
        
        
    def calculation_required(self,atoms, quantities):
        """ Check if a calculation is required. 
        
        Check if the quantities in the quantities list have already been calculated 
        for the atomic configuration atoms. The quantities can be one or more of: 
        'energy', 'forces', and 'stress'.
        """
        raise NotImplementedError


    def set_atoms(self,atoms):
        """ Initialize the calculator for given atomic system. """
        if self.init==False:
            self._initialize(atoms)    
    
    
    def get_occupation_numbers(self,kpt=0,spin=0):
        raise NotImplementedError
        
        
    def get_number_of_bands(self):
        raise NotImplementedError  
               
Hotbit=Calculator    

        
if __name__=='__main__':
    from ase import *
    from box import Atoms
    #h2o=Atoms(symbols='C3H2N2H',positions=[(0,0,0),(2,2,2),(3,3,3),(3,4,5),(1,2,3),(1,2,4),(4,3,2),(4,3,6)],cell=(10,10,10))
    calc=Calculator(SCC=True,txt='test.cal',Anderson_memory=0)
    #h2o.set_calculator(calc)
    #h2o.get_potential_energy()
    #calc.st.solve()
    
    #t=mix.Timer()
    #large=Atoms(symbols='100C',positions=[(i*1.42,0,0) for i in range(100)], cell=(1000,20,20))
    
    
    #large.set_calculator(calc)
    ##t()
    ##for i in range(10):
    #large.get_potential_energy()
    ##t()
    ##for i in range(10):
    
    #calc.st.solve()
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
    calc.st.solve()
    
    
 

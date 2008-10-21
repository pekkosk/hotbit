# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

"""
    ASE-calculator interface for HOTBIT.
    (Hybrid Open-Source Tight-Binding Tool)
    
"""
import numpy as nu
from ase.units import Bohr, Hartree
from box.timing import Timer
from elements import Elements
from interactions import Interactions
from electrostatics import Electrostatics
from repulsion import Repulsion
from states import States
from hotbit.output import Output
import box.mix as mix

         
    
class Calculator(Output):
    """
    ASE-calculator frontend for HOTBIT calculations.
    """
    def __init__(self,parameters=None,
                      elements=None,
                      tables=None,
                      verbose=True,
                      charge=0.0,
                      width=0.02,
                      SCC=True,
                      convergence=1E-3,
                      mixing_constant=0.2,
                      Anderson_memory=3,
                      maxiter=50,
                      gamma_cut=None,
                      txt=None,
                      verbose_SCC=False):
        """ 
        Initialize calculator. 
        
        Parameters:
        -----------
        parameters:       The directory for parametrizations. If parameters==None, use
                          HOTBIT_PARAMETERS environment variable. Parametrizations given
                          by 'elements' and 'tables' keywords override parametrizations
                          in this directory.
         
        elements:         Dictionary for elements used, e.g. {'H':'H_custom.elm','C':'/../C.elm'}
                          Items can also be elements directly: {'H':H} (with H being type Element)
                          * If elements==None, use element info from default directory.
                          * If elements['others']=='default', use default parameters for all other
                            elements than the ones specified. E.g. {'H':'H.elm','others':'default'}
                            (otherwise all elements present have to be specified excplicitly).
                            
        tables:           Dictionary for interactions, e.g. {'CH':'C_H.par','CC':'C_C.par'}
                          * If elements==None, use default interactions.
                          * If elements['others']='default', use default parameters for all other
                            interactions than the ones specified. 
                            E.g. {'CH':'C_H.par','others':'default'}
        
        charge            total electric charge for system (-1 means an additional electron)
        width             width of Fermi occupation (in eV)
        SCC               Self-Consistent Charge calculation
        convergence       convergence criterion (Mulliken charge variation/atom)
        mixing_constant   density mixing for Anderson mixer
        Anderson_memory   Memory for Anderson mixing
        maxiter           Maximum number of self-consistent iterations (for SCC)
        gamma_cut         Range for Coulomb interaction
        txt               Filename for log-file (stdout = None)
        verbose_SCC       Increase verbosity for SCC iterations.
        """       
        from copy import copy
        import os
        
        if gamma_cut!=None: gamma_cut=gamma_cut/Bohr                          
        self.args={ 'parameters':parameters,
                    'elements':elements,
                    'tables':tables,
                    'verbose':verbose,
                    'charge':charge,
                    'width':width/Hartree,
                    'SCC':SCC,
                    'convergence':convergence,
                    'mixing_constant':mixing_constant,
                    'Anderson_memory':Anderson_memory,
                    'maxiter':maxiter,
                    'gamma_cut':gamma_cut,
                    'txt':txt,
                    'verbose_SCC':verbose_SCC}                    
            
        if parameters!=None:
            os.environ.data['HOTBIT_PARAMETERS']=parameters
                        
        self.init=False
        self.element_files=elements
        self.table_files=tables    
        self.verbose=verbose     
        self.set_enabled=True
        self.notes=[]
        self.set_text(self.args['txt'])
        self.timer=Timer('Hotbit',txt=self.get_output())
       
       
    def set(self,key,value):       
        if self.init==True or key not in ['charge']:
            raise AssertionError('Parameters cannot be set after initialization.')
        self.__dict__[key]=value
        
        
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
        from os import uname, popen
        from os.path import abspath, curdir
        from os import environ
        
        revision=popen('svnversion %s' %environ.get('HOTBIT_DIR') ).readline()
        self.version='0.1 (svn=%s)' %revision[:-1]
        print>>self.txt,  '\n\n\n\n\n'
        print>>self.txt,  ' _           _    _     _ _'
        print>>self.txt,  '| |__   ___ | |_ | |__ |_| |_'
        print>>self.txt,  '|  _ \ / _ \|  _||  _ \| |  _|'
        print>>self.txt,  '| | | | ( ) | |_ | ( ) | | |_'
        print>>self.txt,  '|_| |_|\___/ \__|\____/|_|\__|  ver.',self.version
        print>>self.txt,  'Distributed under GNU GPL; see %s' %environ.get('HOTBIT_DIR')+'/LICENSE'
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
       
       
    def get(self,arg=None):
        if arg==None:
            return self.args
        else:
            return self.args[arg]
        
        
    def solve_ground_state(self,atoms):
        """ If atoms moved, solve electronic structure. """
        if not self.init:
            self._initialize(atoms)
                    
        #print 'required?',self.el.calculation_required(atoms,'ground state')
        #print atoms.get_positions()[1]
        #print self.el.atoms.get_positions()[1]
        if self.el.calculation_required(atoms,'ground state'):
            self.el.set_atoms(atoms)
            self.st.solve()
        else:
            pass
                  
                  
    def _initialize(self,atoms):
        """ Initialization of hotbit. """  
        self.timer.start('initialization')
        self.init=True
        self.el=Elements(self,atoms,self.timer,self.element_files,charge=self.args['charge'])
        self.ia=Interactions(self,self.timer,self.el,self.table_files)
        self.es=Electrostatics(self,self.timer)
        self.st=States(self,self.timer,self.el,self.ia)
        self.rep=Repulsion(self,self.timer,self.el,self.ia)
        self.set_enabled=False
        self.pbc=atoms.get_pbc()
        self.el.update_geometry()
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
    
    
    #def calculation_required(self,atoms,quantities):
        #""" Check if a calculation is required. 
        
        #Check if the quantities in the quantities list have already been calculated 
        #for the atomic configuration atoms. The quantities can be one or more of: 
        #'ground state', 'energy', 'forces', and 'stress'.
        #"""
        #return self.el.calculation_required(atoms,quantities)
        
    
    # some not implemented ASE-assumed methods
    def get_fermi_level(self):
        raise NotImplementedError
        

    def set_atoms(self,atoms):
        """ Initialize the calculator for given atomic system. """
        if self.init==True and atoms.get_chemical_symbols()!=self.el.atoms.get_chemical_symbols():
            atoms2=Atoms(atoms)
            raise RuntimeError('Calculator initialized for %s. Create new calculator for %s.'
                               %(self.el.atoms.get_name(),atoms2.get_name() ))
        else:                               
            self._initialize(atoms)    
    
    
    def get_occupation_numbers(self,kpt=0,spin=0):
        raise NotImplementedError
        
        
    def get_number_of_bands(self):
        raise NotImplementedError  
               
Hotbit=Calculator    

     
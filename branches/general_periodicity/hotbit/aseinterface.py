# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

"""
    ASE-calculator interface for HOTBIT.
    (Hybrid Open-Source Tight-Binding Tool)

"""
import numpy as nu
from ase.units import Bohr, Hartree
from ase import Atoms
from box.timing import Timer
from elements import Elements
from interactions import Interactions
from environment import Environment
from repulsion import Repulsion
from states import States
from hotbit.output import Output
import box.mix as mix
from time import time

from atoms import WedgeAtoms



class Calculator(Output):
    """
    ASE-calculator frontend for HOTBIT calculations.
    """
    def __init__(self,parameters=None,
                      elements=None,
                      tables=None,
                      verbose=True,
                      charge=0.0,
                      SCC=True,
                      kpts=(1,1,1),
                      maxiter=50,
                      gamma_cut=None,
                      txt=None,
                      verbose_SCC=False,
                      width=0.02,
                      mixer=None,
                      filename=None):
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

        mixer             Density mixer. Can be a name of the mixer (mixer='Anderson') or a dictionary of parameters (mixer={'name':'Anderson','mixing_constant':0.3, 'memory':5}).
        charge            total electric charge for system (-1 means an additional electron)
        width             width of Fermi occupation (in eV)
        SCC               Self-Consistent Charge calculation
        kpts              number of k-points wrt different symmetries
        maxiter           Maximum number of self-consistent iterations (for SCC)
        gamma_cut         Range for Coulomb interaction
        txt               Filename for log-file (stdout = None)
        verbose_SCC       Increase verbosity for SCC iterations.
        """
        from copy import copy
        import os
        
        if SCC: raise NotImplementedError

        if gamma_cut!=None: gamma_cut=gamma_cut/Bohr
        self.__dict__={ 'parameters':parameters,
                        'elements':elements,
                        'tables':tables,
                        'verbose':verbose,
                        'charge':charge,
                        'width':width/Hartree,
                        'SCC':SCC,
                        'kpts':kpts,
                        'maxiter':maxiter,
                        'gamma_cut':gamma_cut,
                        'txt':txt,
                        'verbose_SCC':verbose_SCC,
                        'mixer':mixer}

        if parameters!=None:
            os.environ.data['HOTBIT_PARAMETERS']=parameters

        self.init=False
        self.notes=[]
        self.set_text(self.txt)
        self.timer=Timer('Hotbit',txt=self.get_output())

        if filename is not None:
            self.load(filename)


    def __copy__(self):
        """
        Returns an uninitialized calculator.
        """
        import sys
        if self.init == True:
            raise NotImplementedError('Calculator has been initialized and it cannot be copied. Please create a new calculator.')

        params_to_copy = ['SCC',
                          'mixer',
                          'width',
                          'tables',
                          'charge',
                          'verbose',
                          'maxiter',
                          'elements',
                          'gamma_cut',
                          'parameters',
                          'verbose_SCC']

        ret = Hotbit()
        ret.__dict__['txt'] = sys.stdout
        for key in params_to_copy:
            ret.__dict__[key] = self.__dict__[key]
        return ret


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
        Output.__del__(self)


    def write(self, filename='restart.hb'):
        """ Write data from calculation to a file. """
        if self.init == False:
            raise RuntimeError('The calculator is not initialized.')
        state = {}
        self.get_data_to_save(state)
        import pickle
        f = open(filename, 'w')
        pickle.dump(state, f)
        f.close()


    def get_data_to_save(self, state):
        """ Gather data from different objects. """
        atoms = {}
        atoms['positions'] = self.el.atoms.get_positions()
        atoms['numbers'] = self.el.atoms.get_atomic_numbers()
        atoms['pbc'] = self.el.atoms.get_pbc()
        atoms['cell'] = self.el.atoms.get_cell()
        state['atoms'] = atoms

        calc = {}
        params = ['parameters','mixer','elements','SCC',
                  'maxiter','tables','gamma_cut','charge','width']
        for key in params:
            calc[key] = self.__dict__[key]
        state['calc'] = calc

        states = {}
        states['prev_dq'] = self.st.prev_dq
        states['count'] = self.st.count
        state['states'] = states


    def load(self, filename):
        """ Load saved state and initialize calculator using that state."""
        import pickle
        f = open(filename)
        state = pickle.load(f)
        f.close()

        atoms = state['atoms']
        pos = atoms['positions']
        num = atoms['numbers']
        pbc = atoms['pbc']
        cell = atoms['cell']
        atoms = Atoms(positions=pos, numbers=num, pbc=pbc, cell=cell)

        calc = state['calc']
        for key in calc:
            self.__dict__[key] = calc[key]
        self._initialize(atoms)

        states = state['states']
        self.st.prev_dq = states['prev_dq']
        self.st.count = states['count']


    def set(self,key,value):
        if self.init==True or key not in ['charge']:
            raise AssertionError('Parameters cannot be set after initialization.')
        self.__dict__[key]=value


    def get_atoms(self):
        """ Return the current atoms object. """
        atoms = self.el.atoms.copy()
        atoms.calc = self
        return atoms


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
        print>>self.txt,  '       Charge=%4.1f' % self.charge
        # TODO: change box and pbc-info here
        print>>self.txt,  '       Box: (Ang)', nu.array(self.el.get_box_lengths())*Bohr
        print>>self.txt,  '       PBC:',self.el.atoms.get_pbc()
        print>>self.txt,  '       Electronic temperature:', self.width*Hartree,'eV'
        mixer = self.st.solver.mixer
        print>>self.txt,  '       Mixer:', mixer.get('name'), 'with memory =', mixer.get('memory'), ', mixing constant =', mixer.get('beta')
        print>>self.txt, self.el.greetings()
        print>>self.txt, self.ia.greetings()
        print>>self.txt, self.rep.greetings()


    def out(self,text):
        print>>self.txt, text
        self.txt.flush()


    def set_text(self,txt):
        """ Set up the output file. """
        if txt is '-':
            self.txt = open('/dev/null','w')
        elif txt is None:
            from sys import stdout
            self.txt=stdout
        else:
            self.txt=open(txt,'a')
        # check if the output of timer must be changed also
        if 'timer' in self.__dict__:
            self.timer.txt = self.txt


    def get(self,arg=None):
        if arg==None:
            return self.__dict__
        else:
            return self.__dict__[arg]


    def solve_ground_state(self,atoms):
        """ If atoms moved, solve electronic structure. """
        
        if not self.init:
            self._initialize(atoms)

        #print 'required?',self.el.calculation_required(atoms,'ground state')
        #print atoms.get_positions()[1]
        #print self.el.atoms.get_positions()[1]
        if self.el.calculation_required(atoms,'ground state'):
            self.el.set_atoms(atoms)
            t0 = time()
            self.st.solve()
            t1 = time()
            if self.verbose:
                print >> self.get_output(), "Solved in %0.2f seconds" % (t1-t0)
        else:
            pass


    def _initialize(self,atoms):
        """ Initialization of hotbit. """
        self.start_timing('initialization')
        self.init=True
        self.el=Elements(self,atoms)
        self.ia=Interactions(self)
        self.st=States(self)
        self.rep=Repulsion(self)
        self.env=Environment(self)
        self.greetings()
        pbc=atoms.get_pbc()
        if self.get('SCC') and nu.any(pbc) and self.get('gamma_cut')==None:
            raise NotImplementedError('SCC not implemented for periodic systems yet (see parameter gamma_cut).')
        if nu.any(pbc) and abs(self.get('charge'))>0.0:
            raise AssertionError('Charged system cannot be periodic.')
        self.flush()
        #self.el.update_geometry()
        self.stop_timing('initialization')


    def get_potential_energy(self,atoms):
        """ Return the potential energy of present system. """
        self.solve_ground_state(atoms)
        self.start_timing('energies')
        ebs=self.get_band_structure_energy(atoms)
        ecoul=self.get_coulomb_energy(atoms)
        erep=self.rep.get_repulsive_energy()
        self.stop_timing('energies')
        return erep+ebs+ecoul
        #return erep
        #return ebs


    def get_forces(self,atoms):
        """ Return the forces of present system. """
        self.solve_ground_state(atoms)
        self.start_timing('forces')
        fbs=self.st.band_structure_forces()
        frep=self.rep.get_repulsive_forces()
        fcoul=self.st.es.gamma_forces() #zero for non-SCC
        self.stop_timing('forces')
        return (fbs+frep+fcoul)*(Hartree/Bohr)
        #return fbs #*(Hartree/Bohr)
#        return frep*(Hartree/Bohr)
    
    
    def get_DOS(self,width=0.1,window=None,npts=201):
        '''
        Return the full density of states, including k-points.
        Zero is the Fermi-level; spin-degeneracy is not counted.
        
        @param width:  Gaussian broadening, in eV
        @param window: energy window 2-tuple, in eV
        @param npts:   number of data points in output
        '''
        self.start_timing('DOS')
        e = self.st.e.copy()*Hartree - self.get_fermi_level()
        flat = e.flatten()
        mn, mx = flat.min(), flat.max()
        if window is not None:
            mn, mx = window
            
        x, y = [],[]
        for a in range(self.el.norb):
            x = nu.concatenate( (x,e[:,a]) )
            y = nu.concatenate( (y,self.st.wk) )
        x=nu.array(x) 
        y=nu.array(y) 
        self.start_timing('broaden')
        dos = mix.broaden(x, y, width=width, N=npts, a=mn, b=mx)
        self.stop_timing('broaden')
        self.stop_timing('DOS')
        return dos


    def get_band_energies(self,kpts):
        '''
        Return band energies for explicitly given list of k-points.
        
        @param kpts: list of k-points; e.g. kpts=[(0,0,0),(pi/2,0,0),(pi,0,0)]
        '''
        

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
        return self.st.es.coulomb_energy()*Hartree


    #def calculation_required(self,atoms,quantities):
        #""" Check if a calculation is required.

        #Check if the quantities in the quantities list have already been calculated
        #for the atomic configuration atoms. The quantities can be one or more of:
        #'ground state', 'energy', 'forces', and 'stress'.
        #"""
        #return self.el.calculation_required(atoms,quantities)


    # some not implemented ASE-assumed methods
    def get_fermi_level(self):
        return self.st.occu.get_mu() * Hartree


#    def set_atoms(self,atoms):
#        """ Initialize the calculator for given atomic system. """
#        if self.init==True and atoms.get_chemical_symbols()!=self.el.atoms.get_chemical_symbols():
#            raise RuntimeError('Calculator initialized for %s. Create new calculator for %s.'
#                               %(self.el.get_name(),mix.parse_name_for_atoms(atoms)))
#        else:
#            self._initialize(atoms)


    def get_occupation_numbers(self,kpt=0,spin=0):
        raise NotImplementedError


    def get_number_of_bands(self):
        raise NotImplementedError


    def start_timing(self, label):
        self.timer.start(label)


    def stop_timing(self, label):
        self.timer.stop(label)


Hotbit=Calculator



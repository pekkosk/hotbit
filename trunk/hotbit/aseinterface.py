# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

"""
    ASE-calculator interface for HOTBIT.
    (Hybrid Open-source Tight-Binding Tool)

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
from grids import Grids
from hotbit.analysis import MullikenAnalysis
from hotbit.analysis import MullikenBondAnalysis
from hotbit.analysis import DensityOfStates
from hotbit.output import Output
import box.mix as mix
from time import time




class Hotbit(Output):
    def __init__(self,parameters=None,
                      elements=None,
                      tables=None,
                      verbose=False,
                      charge=0.0,
                      SCC=True,
                      kpts=(1,1,1),
                      physical_k_points=True,
                      maxiter=50,
                      gamma_cut=None,
                      txt=None,
                      verbose_SCC=False,
                      width=0.02,
                      mixer=None,
                      coulomb_solver=None,
                      filename=None):
        """
        Hotbit -- density-functional tight-binding calculator
                  for atomic simulation environment (ASE).
        
        

        Parameters:
        -----------
        parameters:       The directory for parametrization files. 
                          * If parameters==None, use HOTBIT_PARAMETERS environment variable. 
                          * Parametrizations given by 'elements' and 'tables' keywords 
                            override parametrizations in this directory.

        elements:         Files for element data (*.elm). 
                          example: {'H':'H_custom.elm','C':'/../C.elm'}
                          * Items can also be elements directly: {'H':H} (H is type Element)
                          * If elements==None, use element info from default directory.
                          * If elements['others']=='default', use default parameters for all other
                            elements than the ones specified. E.g. {'H':'H.elm','others':'default'}
                            (otherwise all elements present have to be specified excplicitly).

        tables:           Files for Slater-Koster tables.
                          example: {'CH':'C_H.par','CC':'C_C.par'}
                          * If elements==None, use default interactions.
                          * If elements['others']='default', use default parameters for all other
                            interactions, e.g. {'CH':'C_H.par','others':'default'}

        mixer:            Density mixer. 
                          example: {'name':'Anderson','mixing_constant':0.3, 'memory':5}.
        charge:           Total charge for system (-1 means an additional electron)
        width:            Width of Fermi occupation (eV)
        SCC:              Self-Consistent Charge calculation
                          * True for SCC-DFTB, False for DFTB
        kpts:             Number of k-points.
                          * For translational symmetry points are along the directions
                            given by the cell vectors.
                          * For general symmetries, you need to look at the info
                            from the container used
        physical_k_points Use physical (realistic) k-points for generally periodic systems.
                          * Ignored with normal translational symmetry
                          * True for physically allowed k-points in periodic symmetries.
        maxiter:          Maximum number of self-consistent iterations 
                          * only for SCC-DFTB
        coulomb_solver:   Select Coulomb solver. Choices are
                          * None:     direct summation, see also gamma_cut
                          * 'ewald':  Ewald summation, only works for bravais
                                      containers and becomes slow for large
                                      systems.
                          * 'me':     Multipole expansion
        gamma_cut:        Range for Coulomb interaction if direct summation is
                          selected (see coulomb_solver)
                          * At the moment Coulomb interaction is 1/r*(1-erf(r/gamma_cut))
                            for faster convergence. This is quick & dirty way to
                            calculate electostatics and will be implemented properly.
        txt:              Filename for log-file.
                          * None: standard output
                          * '-': throw output to trash (/null) 
        verbose_SCC:      Increase verbosity in SCC iterations.
        """
        from copy import copy
        import os

        if gamma_cut!=None: gamma_cut=gamma_cut/Bohr
        self.__dict__={ 'parameters':parameters,
                        'elements':elements,
                        'tables':tables,
                        'verbose':verbose,
                        'charge':charge,
                        'width':width/Hartree,
                        'SCC':SCC,
                        'kpts':kpts,
                        'physical_k_points':physical_k_points,
                        'maxiter':maxiter,
                        'gamma_cut':gamma_cut,
                        'txt':txt,
                        'verbose_SCC':verbose_SCC,
                        'mixer':mixer,
                        'coulomb_solver':coulomb_solver}

        if parameters!=None:
            os.environ.data['HOTBIT_PARAMETERS']=parameters

        self.init=False
        self.notes=[]
        #self.set_text(self.txt)
        #self.timer=Timer('Hotbit',txt=self.get_output())

        if filename is not None:
            self.load(filename)


    def __copy__(self):
        """
        Returns an uninitialized calculator.
        """
        import sys
        from os.path import sameopenfile
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
                          'txt',
                          'verbose_SCC',
                          'kpts',
                          'coulomb_solver']

        ret = Hotbit()
        # TODO: if output file already opened (unless /null or stdout)
        # open a file with another name
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
        if self.init:
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
        if key == 'txt':
            self.set_text(value)
        elif self.init==True or key not in ['charge']:
            raise AssertionError('Parameters cannot be set after initialization.')
        else:
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
        print>>self.txt,  '       Container', self.el.container_info()
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
        if txt=='-' or txt=='null':
            self.txt = open('/dev/null','w')
        elif hasattr(txt, 'write'):
            self.txt = txt
        elif txt is None:
            from sys import stdout
            self.txt=stdout
        else:
            self.txt=open(txt,'a')
        # check if the output of timer must be changed also
        if 'timer' in self.__dict__:
            self.timer.txt = self.get_output()


    def get(self,arg=None):
        if arg==None:
            return self.__dict__
        else:
            return self.__dict__[arg]


    def solve_ground_state(self,atoms):
        """ If atoms moved, solve electronic structure. """
        if not self.init:
            self._initialize(atoms)
        if self.el.calculation_required(atoms,'ground state'):
            self.el.update_geometry(atoms)
            t0 = time()
            self.st.solve()
            self.el.set_solved('ground state')
            t1 = time()
            self.flags['Mulliken'] = False
            self.flags['DOS'] = False
            self.flags['bonds'] = False
            if self.verbose:
                print >> self.get_output(), "Solved in %0.2f seconds" % (t1-t0)
        else:
            pass


    def _initialize(self,atoms):
        """ Initialization of hotbit. """
        self.init=True
        self.set_text(self.txt)
        self.timer=Timer('Hotbit',txt=self.get_output())
        self.start_timing('initialization')
        self.el=Elements(self,atoms)
        self.ia=Interactions(self)
        self.st=States(self)
        self.rep=Repulsion(self)
        self.env=Environment(self)
        self.gd=Grids(self)
        pbc=atoms.get_pbc()
        # FIXME: gamma_cut -stuff
        #if self.get('SCC') and nu.any(pbc) and self.get('gamma_cut')==None:
        #    raise NotImplementedError('SCC not implemented for periodic systems yet (see parameter gamma_cut).')
        if nu.any(pbc) and abs(self.get('charge'))>0.0:
            raise AssertionError('Charged system cannot be periodic.')
        self.flush()
        self.el.set_atoms(atoms)
        self.greetings()
        self.flags = {}
        self.flags['Mulliken'] = False
        self.flags['DOS'] = False
        self.flags['bonds'] = False
        self.stop_timing('initialization')
        
        
    def calculation_required(self,atoms,quantities):
        """ Check if a calculation is required.

        Check if the quantities in the quantities list have already been calculated
        for the atomic configuration atoms. The quantities can be one or more of:
        'ground state', 'energy', 'forces', 'magmoms', and 'stress'.
        """
        return self.el.calculation_required(atoms,quantities)


    def get_potential_energy(self,atoms):
        """ Return the potential energy of present system. """
        if self.calculation_required(atoms,['energy']):
            self.solve_ground_state(atoms)
            self.start_timing('energy')
            ebs=self.get_band_structure_energy(atoms)
            ecoul=self.get_coulomb_energy(atoms)
            erep=self.rep.get_repulsive_energy()
            self.epot = ebs + ecoul + erep - self.el.efree*Hartree
            self.stop_timing('energy')
            self.el.set_solved('energy')
        return self.epot.copy()


    def get_forces(self,atoms):
        """ 
        Return forces (in eV/Angstrom)
        
        Ftot = F(band structure) + F(coulomb) + F(repulsion).
        """
        if self.calculation_required(atoms,['forces']):
            self.solve_ground_state(atoms)
            self.start_timing('forces')
            fbs=self.st.get_band_structure_forces()
            frep=self.rep.get_repulsive_forces()
            fcoul=self.st.es.gamma_forces() #zero for non-SCC
            self.stop_timing('forces')
            self.f = (fbs+frep+fcoul)*(Hartree/Bohr)
            self.el.set_solved('forces')
        return self.f.copy()
    
    


    def get_band_energies(self,kpts):
        '''
        Return band energies for explicitly given list of k-points.
        
        @param kpts: list of k-points; e.g. kpts=[(0,0,0),(pi/2,0,0),(pi,0,0)]
        '''
        raise NotImplementerError
        

    def get_stress(self,atoms):
        self.solve_ground_state(atoms)
        return None


    def get_charge(self):
        """ Return system's total charge. """
        return self.get('charge')


    def get_eigenvalues(self,atoms):
        self.solve_ground_state(atoms)
        return self.st.get_eigenvalues()*Hartree


    def get_occupations(self):
        #self.solve_ground_state(atoms)
        return self.st.get_occupations()


    def get_band_structure_energy(self,atoms):
        if self.calculation_required(atoms, ['ebs']):
            self.solve_ground_state(atoms)
            self.ebs = self.st.band_structure_energy()*Hartree
            self.el.set_solved('ebs')
        return self.ebs 


    def get_coulomb_energy(self,atoms):
        if self.calculation_required(atoms,['ecoul']):
            self.solve_ground_state(atoms)
            self.ecoul = self.st.es.coulomb_energy()*Hartree 
            self.st
        return self.ecoul

    def get_grid_density(self,spacing=0.2,pad=False):
        """
        Return the valence electron density
        
        Note that putting wave functions on grid may be very slow; 
        use sparse grids and small cells.
        
        @param spacing: grid spacing. Uses the ASE cell.
        @param pad: if True, the grid spans the whole ASE cell
                    (otherwise one grid point short)
        """
        self.gd.make_grid(spacing=spacing,pad=pad)
        return self.gd.get_density()
    
    
    def get_grid_wf(self,i,k=0,spacing=0.2,pad=False):
        """ 
        Return state i with given k-point on grid. 
        """
        self.gd.make_grid(spacing=spacing,pad=pad)
        return self.gd.get_wf(i)
    

    # some not implemented ASE-assumed methods
    def get_fermi_level(self):
        """
        Return the Fermi-energy (chemical potential) in eV.
        """
        return self.st.occu.get_mu() * Hartree


    def set_atoms(self,atoms):
        """ Initialize the calculator for given atomic system. """
        if self.init==True and atoms.get_chemical_symbols()!=self.el.atoms.get_chemical_symbols():
            raise RuntimeError('Calculator initialized for %s. Create new calculator for %s.'
                               %(self.el.get_name(),mix.parse_name_for_atoms(atoms)))
        else:
            self._initialize(atoms)


    def get_occupation_numbers(self,kpt=0):
        """ Return occupation numbers for given k-point index. """
        return self.st.f[kpt].copy()
    

    def get_number_of_bands(self):
        """ Return the total number of orbitals. """
        return self.st.norb


    def start_timing(self, label):
        self.timer.start(label)


    def stop_timing(self, label):
        self.timer.stop(label)

    #
    # Mulliken population analysis tools
    #
    def _init_mulliken(self):
        """ Initialize Mulliken analysis. """ 
        if self.flags['Mulliken']==False:
            self.MA = MullikenAnalysis(self)
            self.flags['Mulliken']=True
        
    def get_dq(self):
        """ Return atoms' Mulliken populations. """
        return self.st.get_dq()
        
    def get_atom_mulliken(self,I):
        """
        Return Mulliken population for atom I.
        
        parameters:
        ===========
        I:        atom index
        """
        self._init_mulliken()
        return self.MA.atom_mulliken(I)
        
        
    def get_basis_mulliken(self,mu):
        """
        Return Mulliken population of given basis state.
        
        parameters:
        ===========
        mu:     orbital index (see Elements' methods for indices)
        """ 
        self._init_mulliken()
        return self.MA.basis_mulliken(mu)
    
    
    def get_atom_wf_mulliken(self,I,k,a,wk=True):
        """
        Return Mulliken population for given atom and wavefunction.
        
        parameters:
        ===========
        I:      atom index
        k:      k-vector index
        a:      eigenstate index
        wk:     embed k-point weight in population
        """
        self._init_mulliken()
        return self.MA.atom_wf_mulliken(I,k,a,wk)
        
    
    def get_atom_wf_all_orbital_mulliken(self,I,k,a):
        """
        Return orbitals' Mulliken populations for given atom and wavefunction.
        
        parameters:
        ===========
        I:      atom index (returned array size = number of orbitals on I)
        k:      k-vector index 
        a:      eigenstate index
        """
        self._init_mulliken()
        return self.MA.atom_wf_all_orbital_mulliken(I,k,a)
        
    
    def get_atom_wf_all_angmom_mulliken(self,I,k,a,wk=True):
        """ 
        Return atom's Mulliken populations for all angmom for given wavefunction.
        
        parameters:
        ===========
        I:        atom index
        k:        k-vector index
        a:        eigenstate index
        wk:       embed k-point weight into population
        
        return: array (length 3) containing s,p and d-populations      
        """
        self._init_mulliken()
        return self.MA.atom_wf_all_angmom_mulliken(I,k,a,wk)

        
    #
    #  Densities of states methods
    #
    def _init_DOS(self):
        """ Initialize Density of states analysis. """
        if self.flags['DOS']==False:
            self.DOS = DensityOfStates(self)
            self.flags['DOS']=True
            
            
    def get_local_density_of_states(self,projected=False,width=0.05,window=None,npts=501):
        """
        Return state density for all atoms as a function of energy.
        
        parameters:
        ===========
        projected: return local density of states projected for 
                   angular momenta 0,1 and 2 (s,p and d) 
        width:     energy broadening (in eV)
        window:    energy window around Fermi-energy; 2-tuple (eV)
        npts:      number of grid points for energy
        
        return:    projected==False:
                        energy grid, ldos[atom,grid]
                   projected==True:
                        energy grid, 
                        ldos[atom, grid],
                        pldos[atom, angmom, grid]
        """
        self._init_DOS()
        return self.DOS.local_density_of_states(projected,width,window,npts)
    
        
    def get_density_of_states(self,broaden=False,width=0.05,window=None,npts=501):
        """
        Return the full density of states.
        
        Sum of states over k-points. Zero is the Fermi-level.
        Spin-degeneracy is NOT counted.
        
        parameters:
        ===========
        broaden:     If True, return broadened DOS in regular grid
                     in given energy window. 
                     If False, return energies of all states, followed
                     by their weights. 
        width:       Gaussian broadening (eV)
        window:      energy window around Fermi-energy; 2-tuple (eV)
        npts:        number of data points in output
        """
        self._init_DOS()
        return self.DOS.density_of_states(broaden,width,window,npts)
    


    # Bonding analysis
    def _init_bonds(self):
        """ Initialize Mulliken bonding analysis. """
        if self.flags['bonds']==False:
            self.bonds = MullikenBondAnalysis(self)
            self.flags['bonds']=True
            
            
    def get_atom_energy(self,I):
        """ 
        Return the energy of atom I (in eV).
        
        parameters:
        ===========
        I:         atom index
        """
        self._init_bonds(self)
        return self.bonds.atom_energy(I)
    
    
    def get_promotion_energy(self, I):
        """ 
        Return atom's promotion energy (in eV). 
        
        E_prom = sum_(mu in I) q_mu*H^0_(mu,mu)(k) - E_free(I)
        
        parameters:
        ===========
        I:         atom index
        """
        self._init_bonds()
        return self.bonds.promotion_energy(I)
    
    
    def get_bond_energy(self):
        raise NotImplementedError
    
    def get_atom_and_bonding_energy(self):
        raise NotImplementedError

    def get_mayer_bond_order(self):
        raise NotImplementedError
    
    def get_ecov(self):
        raise NotImplementedError
    
    def get_atom_pair_ecov(self):
        """ ecov_I_J"""
        raise NotImplementedError
    
    def get_orbital_pair_ecov(self):
        raise NotImplementedError
    
    def get_angmom_ecov(self):
        raise NotImplementedError
    
 
        

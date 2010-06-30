# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

"""
    ASE-calculator interface for HOTBIT.
    (Hybrid Open-source Tight-Binding Tool)

"""
import os
import glob
import sys
import numpy as np
from auxil import k_to_kappa_points
from ase.units import Bohr, Hartree
from ase import Atoms
from box.timing import Timer
from elements import Elements
from interactions import Interactions
from environment import Environment
from pairpotential import PairPotential
from repulsion import Repulsion
from states import States
from grids import Grids
from hotbit.version import hotbit_version
from hotbit.analysis import MullikenAnalysis
from hotbit.analysis import MullikenBondAnalysis
from hotbit.analysis import DensityOfStates
from hotbit.output import Output
from hotbit.vdw import setup_vdw
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
                      rs='kappa',
                      physical_k=True,
                      maxiter=50,
                      gamma_cut=None,
                      txt=None,
                      verbose_SCC=False,
                      width=0.02,
                      mixer=None,
                      coulomb_solver=None,
                      charge_density='Gaussian',
                      vdw=False,
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
                          * If extension '.elm' is omitted, it is assumed. 
                          * Items can also be elements directly: {'H':H} (H is type Element)
                          * If elements==None, use element info from default directory.
                          * If elements['rest']=='default', use default parameters for all other
                            elements than the ones specified. E.g. {'H':'H.elm','rest':'default'}
                            (otherwise all elements present have to be specified explicitly).

        tables:           Files for Slater-Koster tables.
                          example: {'CH':'C_H.par','CC':'C_C.par'}
                          * If extension '.par' is omitted, it is assumed.
                          * If tables==None, use default interactions.
                          * If tables['rest']='default', use default parameters for all other
                            interactions, e.g. {'CH':'C_H.par','rest':'default'}
                          * If tables['AB']==None, ignore interactions for A and B 
                            (both chemical and repulsive)

        mixer:            Density mixer. 
                          example: {'name':'Anderson','mixing_constant':0.2, 'memory':5}.
        charge:           Total charge for system (-1 means an additional electron)
        width:            Width of Fermi occupation (eV)
        SCC:              Self-Consistent Charge calculation
                          * True for SCC-DFTB, False for DFTB
        kpts:             Number of k-points.
                          * For translational symmetry points are along the directions
                            given by the cell vectors.
                          * For general symmetries, you need to look at the info
                            from the container used
        rs:               * 'kappa': use kappa-points 
                          * 'k': use normal k-points. Only for Bravais lattices.
        physical_k        Use physical (realistic) k-points for generally periodic systems.
                          * Ignored with normal translational symmetry
                          * True for physically allowed k-points in periodic symmetries.
        maxiter:          Maximum number of self-consistent iterations 
                          * only for SCC-DFTB
        coulomb_solver:   The Coulomb solver object. If None, a DirectCoulomb
                          object will the automatically instantiated.
                          * only for SCC-DFTB
        charge_density:   Shape of the excess charge on each atom. Possibilities
                          are:
                          * 'Gaussian': Use atom centered Gaussians. This is the
                            default.
                          * 'Slater': Slater-type exponentials as used in the
                            original SCC-DFTB scheme.
                          * only for SCC-DFTB
        gamma_cut:        Range for Coulomb interaction if direct summation is
                          selected (coulomb_solver = None).
                          * only for SCC-DFTB
        vdw               Include van der Waals interactions
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
                        'rs':rs,
                        'physical_k':physical_k,
                        'maxiter':maxiter,
                        'gamma_cut':gamma_cut,
                        'vdw':vdw,
                        'txt':txt,
                        'verbose_SCC':verbose_SCC,
                        'mixer':mixer,
                        'coulomb_solver':coulomb_solver,
                        'charge_density':charge_density}

        if parameters!=None:
            os.environ.data['HOTBIT_PARAMETERS']=parameters

        self.init=False
        self.notes=[]
        self.dry_run = '--dry-run' in sys.argv
        #self.set_text(self.txt)
        #self.timer=Timer('Hotbit',txt=self.get_output())

        if filename is not None:
            self.load(filename)

    def copy(self):
        return self.__copy__()

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
                          'rs',
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
        
        
    def write_electronic_data(self,filename,keys=None):
        """
        Write key electronic data into a file with *general* format.
        
        Hotbit is not needed to analyze the resulting data file.
        The data will be in a dictionary with the following items:
        
        N          the number of atoms
        norb       the number of orbitals
        nelectrons the number of electrons
        charge     system charge
        epot       potential energy
        ebs        band structure energy
        ecoul      coulomb energy
        erep       repulsive energy
        forces     atomic forces
        symbols    element symbols
        e          single-particle energies
        occ        occupations
        nk         number of k-points
        k          k-point vectors
        wk         k-point weights
        dq         excess Mulliken populations
        gap        energy gap
        gap_prob   certainty of the gap determination above
        dose       energies for density of states (all states over k-points as well)
                   0 = Fermi-level
        dos        density of states (including k-point weights)
        
        Access to data, simply:
        
        data = numpy.load(filename)
        print data['epot'] 
        
        parameters:
        -----------
        filename:     output file name
        keys:         list of items (key names) to save. 
                      If None, save all.
        """ 
        data = {}
        data['N'] = self.el.N
        data['norb'] = self.st.norb
        data['charge'] = self.get('charge')
        data['nelectrons'] = self.el.get_number_of_electrons()
        data['erep'] = self.rep.get_repulsive_energy()
        data['ecoul'] = self.get_coulomb_energy(self.el.atoms)
        data['ebs'] = self.get_band_structure_energy(self.el.atoms)
        data['epot'] = self.get_potential_energy(self.el.atoms)
        data['forces'] = self.get_forces(self.el.atoms)
        data['symbols'] = self.el.symbols
        data['e'] = self.st.e
        data['occ'] = self.st.f
        data['nk'] = self.st.nk
        data['k'] = self.st.k
        data['wk'] = self.st.wk
        data['dq'] = self.st.mulliken()
        data['gap'], data['gap_prob'] = self.get_energy_gap()
        data['dose'], data['dos'] = self.get_density_of_states(False)
        
        for key in data.keys():
            if keys!=None and key not in keys:
                del data[key]
        import pickle  
        f = open(filename, 'w')
        pickle.dump(data,f)
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
        params = ['parameters','mixer','elements','SCC','rs',
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

        self.version=hotbit_version
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
        print>>self.txt,  'Symmetry operations (if any):'
        rs = self.get('rs')
        kpts = self.get('kpts')
        if not isinstance(kpts,tuple):
            kpts = (NaN,NaN,NaN)
        M = self.el.get_number_of_transformations()
        for i in range(3):
            print>>self.txt, '       %i: pbc=' %i, self.el.atoms.get_pbc()[i],
            print>>self.txt, ', %s-points=%i, M=%.f' %(rs,kpts[i],M[i])
        print>>self.txt,  'Electronic temperature:', self.width*Hartree,'eV'
        mixer = self.st.solver.mixer
        print>>self.txt,  'Mixer:', mixer.get('name'), 'with memory =', mixer.get('memory'), ', mixing constant =', mixer.get('beta')
        print>>self.txt, self.el.greetings()
        print>>self.txt, self.ia.greetings()
        print>>self.txt, self.rep.greetings()
        if self.pp.exists():
            print>>self.txt, self.pp.greetings()


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
        
        
    def memory_estimate(self):
        """
        Return an estimate for memory consumption in GB.
        
        If script run with --dry-run, print estimate and exit.
        """
        if self.st.nk>1:
            number = 16. #complex
        else:
            number = 8 #real
        M = self.st.nk*self.st.norb**2*number
        #     H   S   dH0   dS    wf  H1  dH   rho rhoe
        mem = M + M + 3*M + 3*M + M + M + 3*M + M + M
        if self.dry_run:
            print 'Memory consumption estimate %.2f GB' %(mem/1E9) 
            raise SystemExit
        else:
            return mem/1E9       


    def solve_ground_state(self,atoms):
        """ If atoms moved, solve electronic structure. """
        if not self.init:
            self._initialize(atoms)
        if self.calculation_required(atoms,'ground state'):
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
            if self.get('SCC'):
                atoms.set_charges(-self.st.get_dq())
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
        self.pp=PairPotential(self)
        if self.get('vdw'):
            setup_vdw(self)
        self.env=Environment(self)
        pbc=atoms.get_pbc()
        # FIXME: gamma_cut -stuff
        #if self.get('SCC') and np.any(pbc) and self.get('gamma_cut')==None:
        #    raise NotImplementedError('SCC not implemented for periodic systems yet (see parameter gamma_cut).')
        if np.any(pbc) and abs(self.get('charge'))>0.0:
            raise AssertionError('Charged system cannot be periodic.')
        self.flush()
        self.el.set_atoms(atoms)
        self.greetings()
        self.flags = {}
        self.flags['Mulliken'] = False
        self.flags['DOS'] = False
        self.flags['bonds'] = False
        self.flags['grid'] = False
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
            epp=self.pp.get_energy()
            self.epot = ebs + ecoul + erep + epp - self.el.efree*Hartree
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
            fpp = self.pp.get_forces()
            self.stop_timing('forces')
            self.f = (fbs+frep+fcoul+fpp)*(Hartree/Bohr)
            self.el.set_solved('forces')
        return self.f.copy()
    

    def get_band_energies(self, kpts=None, shift=True, rs='kappa'):
        '''
        Return band energies for explicitly given list of k-points.
        
        parameters:
        ===========
        kpts:      list of k-points; e.g. kpts=[(0,0,0),(pi/2,0,0),(pi,0,0)]
                   k- or kappa-points, depending on parameter rs.
        shift:     shift zero to the Fermi-level
        rs:        use 'kappa'- or 'k'-points in reciprocal space
        '''
        if kpts==None:
            e = self.st.e * Hartree
        else:
            if rs=='k':
                klist = k_to_kappa_points(kpts,self.el.atoms)
            elif rs=='kappa':
                klist = kpts
            e = self.st.get_band_energies(klist)*Hartree
        
        if shift:
            return e-self.get_fermi_level()
        else:
            return e 


    def get_stress(self,atoms):
        self.solve_ground_state(atoms)
        return None


    def get_charge(self):
        """ Return system's total charge. """
        return self.get('charge')


    def get_eigenvalues(self):
        return self.st.get_eigenvalues()*Hartree
    
    
    def get_energy_gap(self):
        """
        Return the energy gap. (in eV)
        
        Gap is the energy difference between the first states
        above and below Fermi-level. Return also the probability
        of having returned the gap; it is the difference
        in the occupations of these states, divided by 2.
        """
        eigs = (self.get_eigenvalues() - self.get_fermi_level()).flatten()
        occ = self.get_occupations().flatten()
        ehi, elo=1E10,-1E10
        for e,f in zip(eigs,occ):
            if elo<e<=0.0:
                elo = e
                flo = f
            elif 0.0<e<ehi:
                ehi = e
                fhi = f
        return ehi-elo, (flo-fhi)/2      



    def get_occupations(self):
        #self.solve_ground_state(atoms)
        return self.st.get_occupations()


    def get_band_structure_energy(self,atoms):
        if self.calculation_required(atoms, ['ebs']):
            self.solve_ground_state(atoms)
            self.ebs = self.st.get_band_structure_energy()*Hartree
            self.el.set_solved('ebs')
        return self.ebs 


    def get_coulomb_energy(self,atoms):
        if self.calculation_required(atoms,['ecoul']):
            self.solve_ground_state(atoms)
            self.ecoul = self.st.es.coulomb_energy()*Hartree 
            self.st
        return self.ecoul


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
    #   grid stuff
    #
    def set_grid(self,h=0.2,cutoff=3.0):
        if self.calculation_required(self.el.atoms,['energy']):
            raise AssertionError('Electronic structure is not solved yet!')
        if self.flags['grid']==False:
            self.gd = Grids(self,h,cutoff)
            self.flags['grid']=True
        
        
    def get_grid_basis_orbital(self,I,otype,k=0,pad=True):
        """
        Return basis orbital on grid.
        
        parameters:
        ===========
        I:     atom index
        otype: orbital type ('s','px','py',...)
        k:     k-point index (basis functions are really the extended
               Bloch functions for periodic systems)
        pad:   padded edges in the array
        """
        if self.flags['grid']==False:
            raise AssertionError('Grid needs to be set first by method "set_grid".')
        return self.gd.get_grid_basis_orbital(I,otype,k,pad)


    def get_grid_wf(self,a,k=0,pad=True):
        """ 
        Return eigenfunction on a grid.
        
        parameters:
        ===========
        a:     state (band) index
        k:     k-vector index
        pad:   padded edges 
        """
        if self.flags['grid']==False:
            raise AssertionError('Grid needs to be set first by method "set_grid".')
        return self.gd.get_grid_wf(a,k,pad)
    
    
    def get_grid_wf_density(self,a,k=0,pad=True):
        """
        Return eigenfunction density.
        
        Density is not normalized; accurate quantitative analysis
        on this density are best avoided.
        
        parameters:
        ===========
        a:     state (band) index
        k:     k-vector index
        pad:   padded edges
        """
        if self.flags['grid']==False:
            raise AssertionError('Grid needs to be set first by method "set_grid".')
        return self.gd.get_grid_wf_density(a,k,pad)
    
    
    def get_grid_density(self,pad=True):
        """ 
        Return electron density on grid.
        
        Do not perform accurate analysis on this density.
        Integrated density differs from the total number of electrons.
        Bader analysis inaccurate.
        
        parameters:
        pad:      padded edges
        """
        if self.flags['grid']==False:
            raise AssertionError('Grid needs to be set first by method "set_grid".')
        return self.gd.get_grid_density(pad)
            

    #
    # Mulliken population analysis tools
    #
    def _init_mulliken(self):
        """ Initialize Mulliken analysis. """ 
        if self.calculation_required(self.el.atoms,['energy']):
            raise AssertionError('Electronic structure is not solved yet!')
        if self.flags['Mulliken']==False:
            self.MA = MullikenAnalysis(self)
            self.flags['Mulliken']=True
        
    def get_dq(self):
        """ Return atoms' excess Mulliken populations.
        
        The total populations subtracted by
        the numbers of valence electrons.
        
        """
        return self.st.get_dq()
        
    def get_atom_mulliken(self,I):
        """
        Return Mulliken population for atom I.
        
        This is the total population, without the number
        of valence electrons subtracted.
        
        parameters:
        ===========
        I:        atom index
        """
        self._init_mulliken()
        return self.MA.get_atom_mulliken(I)
        
        
    def get_basis_mulliken(self,mu):
        """
        Return Mulliken population of given basis state.
        
        parameters:
        ===========
        mu:     orbital index (see Elements' methods for indices)
        """ 
        self._init_mulliken()
        return self.MA.get_basis_mulliken(mu)
    
    
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
        return self.MA.get_atom_wf_mulliken(I,k,a,wk)
        
    
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
        return self.MA.get_atom_wf_all_orbital_mulliken(I,k,a)
        
    
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
        return self.MA.get_atom_wf_all_angmom_mulliken(I,k,a,wk)

        
    #
    #  Densities of states methods
    #
    def _init_DOS(self):
        """ Initialize Density of states analysis. """
        if self.calculation_required(self.el.atoms,['energy']):
            raise AssertionError('Electronic structure is not solved yet!')
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
        return self.DOS.get_local_density_of_states(projected,width,window,npts)
    
        
    def get_density_of_states(self,broaden=False,projected=False,occu=False,width=0.05,window=None,npts=501):
        """
        Return the full density of states.
        
        Sum of states over k-points. Zero is the Fermi-level.
        Spin-degeneracy is NOT counted.
        
        parameters:
        ===========
        broaden:     * If True, return broadened DOS in regular grid
                       in given energy window. 
                     * If False, return energies of all states, followed
                       by their k-point weights.
        projected:   project DOS for angular momenta 
        occu:        for not broadened case, return also state occupations
        width:       Gaussian broadening (eV)
        window:      energy window around Fermi-energy; 2-tuple (eV)
        npts:        number of data points in output
        
        return:      * if projected: e[:],dos[:],pdos[l,:] (angmom l=0,1,2)
                     * if not projected: e[:],dos[:]
                       * if broaden: e[:] is on regular grid, otherwise e[:] are
                         eigenvalues and dos[...] corresponding weights
                     * if occu: e[:],dos[:],occu[:] 
                          
        """
        self._init_DOS()
        return self.DOS.get_density_of_states(broaden,projected,occu,width,window,npts)
    


    # Bonding analysis
    def _init_bonds(self):
        """ Initialize Mulliken bonding analysis. """
        if self.calculation_required(self.el.atoms,['energy']):
            raise AssertionError('Electronic structure is not solved yet!')
        if self.flags['bonds']==False:
            self.bonds = MullikenBondAnalysis(self)
            self.flags['bonds']=True
            
            
    def get_atom_energy(self,I=None):
        """ 
        Return the energy of atom I (in eV).
        
        Warning: bonding & atom energy analysis less clear for
        systems where orbitals overlap with own periodic images.
        
        parameters:
        ===========
        I:         atom index. If None, return all atoms' energies
                   as an array.
        """
        self._init_bonds()
        return self.bonds.get_atom_energy(I)



    def get_mayer_bond_order(self,i,j):
        """
        Return Mayer bond-order between two atoms.
        
        Warning: bonding & atom energy analysis less clear for
        systems where orbitals overlap with own periodic images.
        
        parameters:
        ===========
        I:        first atom index
        J:        second atom index
        """
        self._init_bonds()
        return self.bonds.get_mayer_bond_order(i,j)


    def get_promotion_energy(self,I=None):
        """ 
        Return atom's promotion energy (in eV). 
        
        Defined as:
            E_prom,I = sum_(mu in I) [q_(mu) - q_(mu)^0] epsilon_mu
        
        parameters:
        ===========
        I:         atom index. If None, return all atoms' energies
                   as an array.
        """
        self._init_bonds()
        return self.bonds.get_promotion_energy(I)


    def get_bond_energy(self,i,j):
        """ 
        Return the absolute bond energy between atoms (in eV). 
        
        Warning: bonding & atom energy analysis less clear for
        systems where orbitals overlap with own periodic images.
        
        parameters:
        ===========
        i,j:     atom indices
        """
        self._init_bonds()
        return self.bonds.get_bond_energy(i,j)
        
    
    def get_atom_and_bond_energy(self,i=None):
        """
        Return given atom's contribution to cohesion.
        
        parameters:
        ===========
        i:    atom index. If None, return all atoms' energies
              as an array.
        """
        self._init_bonds()
        return self.bonds.get_atom_and_bond_energy(i)
        
        
    def get_covalent_energy(self,mode='default',i=None,j=None,width=None,window=None,npts=501):
        """
        Return covalent bond energies in different modes. (eV)
        
        ecov is described in 
        Bornsen, Meyer, Grotheer, Fahnle, J. Phys.:Condens. Matter 11, L287 (1999) and
        Koskinen, Makinen Comput. Mat. Sci. 47, 237 (2009)
        
        
        
        parameters:
        ===========
        mode:    'default' total covalent energy
                 'orbitals' covalent energy for orbital pairs
                 'atoms' covalent energy for atom pairs
                 'angmom' covalent energy for angular momentum components
        i,j:     atom or orbital indices, or angular momentum pairs
        width:   * energy broadening (in eV) for ecov
                 * if None, return energy eigenvalues and corresponding 
                   covalent energies in arrays, directly
        window:  energy window (in eV wrt Fermi-level) for broadened ecov
        npts:    number of points in energy grid (only with broadening) 
    
        return:
        =======
        x,y:     * if width==None, x is list of energy eigenvalues (including k-points)
                   and y covalent energies of those eigenstates
                 * if width!=None, x is energy grid for ecov.
                 * energies (both energy grid and ecov) are in eV.
         
        Note: energies are always shifted so that Fermi-level is at zero. 
              Occupations are not otherwise take into account (while k-point weights are)
        """
        self._init_bonds()
        return self.bonds.get_covalent_energy(mode,i,j,width,window,npts)
    
 
    def add_pair_potential(self,i,j,v,eVA=True):
        """
        Add pair interaction potential function for elements or atoms
        
        parameters:
        ===========
        i,j:    * atom indices, if integers (0,1,2,...)
                * elements, if strings ('C','H',...)
        v:      Pair potential function. 
                Only one potential per element and atom pair allowed. 
                Syntax:  v(r,der=0), v(r=None) returning the
                interaction range in Bohr or Angstrom.
        eVA:    True for v in eV and Angstrom
                False for v in Hartree and Bohr
        """
        self.pp.add_pair_potential(i,j,v,eVA)
        


### Helper functions

def database_from_path(path):
    if path is None:
        path = '.'

    fns = glob.glob('%s/*.elm' % path)

    elements = { }
    tables = { }

    if len(fns) > 0:
        for fn in fns:
            i0 = fn.rfind('/')
            i1 = fn.rfind('.')
            el1 = fn[i0+1:i1]

            elements[el1] = fn

            for fn2 in fns:
                i0 = fn.rfind('/')
                i1 = fn.rfind('.')
                el2 = fn2[i0+1:i1]

                if os.path.exists('%s/%s_%s.par' % ( path, el1, el2 )):
                    tables['%s%s' % ( el1, el2 )] = \
                        '%s/%s_%s.par' % ( path, el1, el2 )
                else:
                    if os.path.exists('%s/%s_%s.par' % ( path, el2, el1 )):
                        tables['%s%s' % ( el1, el2 )] = \
                            '%s/%s_%s.par' % ( path, el2, el1 )

    else:
        fns = glob.glob('%s/*-*.skf' % path)

        if len(fns) > 0:
            for fn in fns:
                i0 = fn.rfind('/')
                i1 = fn.rfind('-')
                i2 = fn.rfind('.')
                el1 = fn[i0+1:i1]
                el2 = fn[i1+1:i2]

                if el1 == el2:
                    elements[el1] = fn
                tables['%s%s' % ( el1, el2 )] = fn
        else:
            raise RuntimeError('No Slater-Koster database found in directory '
                               '%s.' % path)

    return { 'elements': elements, 'tables': tables }

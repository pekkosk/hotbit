from element import Element
from ase import Atoms as ase_Atoms
from hotbit.atoms import Atoms
import numpy as nu
import hotbit.auxil as aux
import box.mix as mix
from numpy.linalg.linalg import norm
from ase.units import Hartree,Bohr
from os import environ,path
from weakref import proxy
from copy import copy, deepcopy
#from fortran.misc import fortran_doublefor



class Elements:
    def __init__(self,calc,atoms,charge=None):
        """
        Read element info from atoms and set up the input files
        for reading element data.

        If atoms is ASE.atoms object, convert it into box.Atoms object
        """
        elements = calc.get('elements')
        if elements!=None:
            elements=elements.copy()

        if not isinstance(atoms, ase_Atoms):
            raise AssertionError('Given atoms object has to be ase.Atoms type.')

        self._update_atoms(atoms)

        self.calc=proxy(calc)
        self.symbols=atoms.get_chemical_symbols()
        self.N=len(atoms)
        self.name = None
        self.positions=None
        if charge == None:
            self.charge = calc.get_charge()

        self.present=[]
        for element in self.symbols:
            if element not in self.present:
                self.present.append(element)

        # default input files if defined. Override them by the custom
        # input files, if specified.
        
        self.files={}
        for key in self.symbols: 
            self.files[key]=None
            
        # set customized files
        current = path.abspath('.')
        default = environ.get('HOTBIT_PARAMETERS')
        if elements!=None:
            for key in elements:
                if key=='rest': continue
                file = elements[key]
                if file[-4:]!='.elm':
                    file+='.elm'
                if not path.isfile(file):
                    raise RuntimeError('Custom element file "%s" for %s not found.' %(file,key))
                else:
                    file = path.abspath(file)
                    self.files[key] = file
        
        # find element data from default place
        if elements==None or elements!=None and 'rest' in elements and elements['rest']=='default':
            for key in self.symbols:
                if self.files[key]!=None: continue
                file = path.join(default,'%s.elm' %key)
                if not path.isfile(file):
                    raise RuntimeError('Default element file "%s" for %s not found.' %(file,key))
                else:
                    self.files[key] = file 
        
        self._elements_initialization()
        self.solved={'ground state':None,'energy':None,'forces':None,'stress':None,'ebs':None,'ecoul':None,'magmoms':None}

    def __len__(self):
        return len(self.symbols)

    def __del__(self):
        pass
    
    def get_N(self):
        """ Return the number of atoms. """
        return self.N    

    def greetings(self):
        """ Return documentation for elements from .elm files. """
        txt='%i atoms, %i states, %.1f electrons (%.1f filled states)\n' %(self.N,self.norb,self.electrons,self.electrons/2)
        txt+=self.range_message+'\n'
        for s in self.present:
            el=self.elements[s]
            comment=el.get_comment()
            if len(comment) > 0:
                if type(self.files[s])==type(''):
                    file=self.files[s]
                    txt+='Element %s in %s\n' %(s,file)
                else:
                    txt+='Element %s (object given)\n' %s
                for line in comment:
                    txt+='    *'+line.lstrip()+'\n'
        if txt=='':
            txt='No comments for elements.'
        return txt    
    
    def set_atoms(self,atoms):
        """ Set the atoms object ready for calculations. """
        self._update_atoms(atoms)
            
        # determine ranges if they go to infinity
        r = self.atoms.get_symmetry_operation_ranges()
        Mlarge = 5 # TODO: chek Mlarge to be large enough
        self.ranges = []
        s = 'Initial n ranges:'
        for i in range(3):
            assert r[i,0]<=r[i,1]
            if not self.atoms.get_pbc()[i]: assert r[i,0]==r[i,1]==0
            if r[i,0]==-nu.Inf:
                assert r[i,1]==nu.Inf
                r[i,:] = [-Mlarge,Mlarge]
            elif r[i,0]<-Mlarge:
                assert r[i,1]>=Mlarge
                r[i,:] = [-Mlarge,Mlarge]

            a,b = int(round(r[i,0])), int(round(r[i,1]))
            s += '[%i,%i] ' %(a,b)
            self.ranges.append( range(a,b+1) )
        self.range_message = s



    def calculation_required(self,atoms,quantities):
        """ Return True if quantities are solved for atoms.

        The quantities can be one or more of: 'ground state', 'energy', 'forces', and 'stress'.
        'ebs' or 'ecoul', 'magmoms'
        """
        if not isinstance(quantities,(list,tuple)):
            quantities=[quantities]

        if type(self.solved['ground state']) == type(None):
            return True

        # check that all quantities have been solved for identical atoms
        for quantity in quantities:
            solved_atoms = self.solved[quantity]
            if type(solved_atoms)==type(None):
                return True                 
            if solved_atoms!=atoms:
                return True
        return False


    def set_solved(self,quantities):
        """ Set quantities solved for current atoms. """
        if not isinstance(quantities,(list,tuple)):
            quantities=[quantities]
        for quantity in quantities:
            self.solved[quantity]=self.atoms.copy()
            

    def _update_atoms(self,atoms):
        """ Update atoms-object, whether it is ase.Atoms or hotbit.Atoms. """
        
        if hasattr(atoms,'container'):
            self.atoms = atoms.copy()
        else:
            bravais = Atoms(atoms=atoms,container='Bravais')
            self.atoms = bravais.copy()


    def update_geometry(self,atoms):
        '''
        Update all properties related to geometry (calculate once/geometry)
        '''
        self.calc.start_timing('geometry')    
        self._update_atoms(atoms) 
        # select the symmetry operations n where atoms still interact chemically
        nmax = len(self.ranges[0])*len(self.ranges[1])*len(self.ranges[2])
        ijn = nu.zeros( (self.N,self.N,nmax),int )
        ijnn = nu.zeros( (self.N,self.N),int )
        
        # add the n=(0,0,0) first (separately)
        self.Rn = [[self.nvector(r=i,ntuple=(0,0,0)) for i in xrange(self.N)]]
        self.Rot = [ self.rotation((0,0,0)) ]
        self.ntuples = [(0,0,0)]
                
        for i in xrange(self.N):
            for j in xrange(self.N):
                dij = nu.linalg.norm( self.Rn[0][i]-self.Rn[0][j] )
                if dij < self.calc.ia.hscut[i,j]:
                    ijn[i,j,ijnn[i,j]] = 0
                    ijnn[i,j] += 1 
        
        
        # calculate the distances from unit cell 0 to ALL other possible; select chemically interacting
        self.calc.start_timing('operations')
        # FIXME!!! This does not consider 'gamma_cut'!
        cut2 = self.calc.ia.hscut**2
        n = 1
        for n1 in self.ranges[0]:
            for n2 in self.ranges[1]:
                for n3 in self.ranges[2]:
                    nt = (n1,n2,n3)
                    if nt==(0,0,0): continue
                    # check that any atom interacts with this unit cell
                    R = nu.array([self.nvector(r=i,ntuple=nt) for i in xrange(self.N)])
                    rn  = nu.array(self.Rn[0])

                    dRt = rn[:, 0].reshape(1,-1)-R[:, 0].reshape(-1,1)
                    dR  = dRt*dRt
                    dRt = rn[:, 1].reshape(1,-1)-R[:, 1].reshape(-1,1)
                    dR += dRt*dRt
                    dRt = rn[:, 2].reshape(1,-1)-R[:, 2].reshape(-1,1)
                    dR += dRt*dRt
                    addn = nu.any(dR <= cut2)
                    
                    if addn:
                        n += 1
                        self.ntuples.append(nt)
                        self.Rn.append(R)
                        self.Rot.append( self.rotation(nt) )
            
        self.ijn = ijn
        self.ijnn = ijnn
        self.Rn = nu.array(self.Rn)
        self.Rot = nu.array(self.Rot)
        self.calc.stop_timing('operations')
        
        self.calc.start_timing('displacements')
        rijn = nu.zeros((len(self.ntuples),self.N,self.N,3))
        dijn = nu.zeros((len(self.ntuples),self.N,self.N))
        # Fixme!!! Think about how to use numpy for speedup!
        for i in range(self.N):
            for j in range(self.N):
                rijn[:,i,j,:] = self.Rn[:,j,:] - self.Rn[0,i,:]
                dijn[:,i,j]   = nu.sqrt( (rijn[:,i,j,:]**2).sum(axis=1) )
        self.rijn = rijn
        self.dijn = dijn
        self.calc.stop_timing('displacements')
                        
        
#===============================================================================
#        # TODO                
#        def check_too_close_distances(self):
#        # FIXME: move this to elements--it's their business; and call in geometry update
#        """ If some element pair doesn't have repulsive potential,
#            check that they are not too close to each other. """
#        for si in self.present:
#            for sj in self.present:
#                d = self.calc.el.distance_of_elements(si,sj,mode='minimum')
#                if d != None and self.kill_radii[si,sj] != None:
#                    if d < self.kill_radii[si,sj]:
#                        raise AssertionError("Atoms with no repulsive potential are too close to each other: %s and %s" % (si, sj))
#===============================================================================

        
        # TODO: calc.ia should also know the smallest allowed distances between elements
        # (maybe because of lacking repulsion or SlaKo tables), this should be checked here!

        self.calc.stop_timing('geometry')
        
    def get_pbc(self):
        return self.atoms.get_pbc()

    def get_transforms(self):
        return self.ntuples

    def get_distances(self):
        # FIXME!!! Update distances first?
        return self.rijn, self.dijn

       
    def rotation_of_axes(self,n):
        '''
        Return the quantization axis rotation matrix for given symmetry operation.
         
        @param n: 3-tuple for transformation 
        '''
        return self.atoms.rotation_of_axes(n)
    
    
    def rotation(self,n):
        '''
        Return the quantization axis rotation matrix for given symmetry operation.
         
        @param n: 3-tuple for transformation 
        '''
        return self.atoms.rotation(n)
        
    
    def nvector(self,r,ntuple=(0,0,0),r0=nu.array([0,0,0]),lst='vec'):
        '''
        Return position vector rn-r0, when r is operated by S(n)r=rn.
        
        Positions should be in atomic units.
        @param r:   position (array) or atom index (integer) which is operated
        @param ntuple:   operate on r with S(n)
        @param r0:  if integer, use atom r0's position as r0
        @param l:   list of properties to return, 'vec'=vector, 'hat'=unit vector, 'norm'=norm
        '''
        if not isinstance(lst,(list,tuple)):
            lst=[lst]
        assert not( nu.all(r0!=nu.array([0,0,0])) and 'tensor' in lst )
        if isinstance(r,int):  
            r=self.atoms.positions[r] 
        else:
            r=r*Bohr
        if isinstance(r0,int): 
            r0=self.atoms.positions[r0] 
        else:
            r0=r0*Bohr
        vec=(self.atoms.transform(r,ntuple)-r0)/Bohr
        
        ret=[]
        for l in lst:
            if l=='vec': 
                ret.append(vec)
            elif l=='hat':
                norm = nu.linalg.norm(vec)
                if norm<1E-6: raise AssertionError('Suspiciously short vector')
                ret.append( vec/norm )
            elif l=='norm': 
                ret.append( nu.linalg.norm(vec) )
            else:
                raise AssertionError('Keyword %s not defined' %l)
        
        if len(ret)==1: 
            return ret[0]
        else: 
            return ret

    def set_cutoffs(self,cutoffs):
        """ Set the maximum interaction cutoffs (dict of ranges for element pair SlaKo tables). """
        self.cut=cutoffs
        self.maxcut=0.0
        for key in self.cut:
            self.maxcut=max(self.maxcut,self.cut[key])


    def get_name(self):
        """ Get the name of the system, e.g. 'benzene'. Default are symbols, 'H2O' """
        if self.name == None:
            self.name = mix.parse_name_for_atoms(self.atoms)
        return self.name
    
    def container_info(self):
        return repr(self.atoms.container)

    def set_name(self, name):
        self.name = name


    def _elements_initialization(self):
        '''
        Initialize element objects, orbital tables etc.
        
        This initialization is done only once for given set of element info.
        Initialization of any geometrical properties is done elsewhere. 
        '''
        self.elements={}
        for symb in self.present:
            if symb not in self.files:
                raise KeyError('Element file for %s was not defined.' %symb)
            if type(self.files[symb])==type(''):
                self.elements[symb]=Element(self.files[symb])
            else:
                self.elements[symb]=self.files[symb]
        self.efree = self.get_free_atoms_energy()

        self.orb=[]
        self.pbc = self.atoms.get_pbc()
        # index: the number of the orbital
        # atomindex: the number of the orbital on the atom
        # atom: the number of atom the orbital is centered on
        self.atomorb=[]
        self.nr_orbitals=[]
        self.first_orbitals=[]
        self.orbital_atoms=[]
        self.norb=0
        ls = ['s','p','d']
        for i,symb in enumerate(self.symbols):
            el=self.elements[symb]
            atomorb=[]
            self.first_orbitals.append(self.norb)
            for k,(ao,e) in enumerate(zip(el.get_orbital_types(),el.get_onsite_energies())):
                self.norb+=1
                Rnl=el.get_Rnl_function(ao)
                self.orbital_atoms.append(i)
                angmom = ls.index(ao[0])
                atomorb.append({'atom':i,'symbol':symb,'orbital':ao,'angmom':angmom,\
                                'index':self.norb-1,'energy':e,'atomindex':k,'Rnl':Rnl})
            self.nr_orbitals.append(len(atomorb))
            self.atomorb.append(atomorb)
            self.orb.extend(atomorb)

        # number of valence electrons on atoms and the total number of electrons
        self.nr_of_valences=nu.array( [self.elements[symb].get_valence_number() for symb in self.symbols] )
        self.electrons=self.get_valences().sum()-self.charge

        # set list of element pairs
        #self.el_pair_list=[]
        #for i,si in enumerate(self.present):
        #    for j,sj in enumerate(self.present):
        #        self.el_pair_list.append([si,sj])

        self.atom_orb_indices=[[orb['index'] for orb in self.atomorb[i]] for i in range(self.N)]

        # array version for fortran manipulations
        self.atom_orb_indices2=nu.zeros((self.N,9),int)-1
        for i,noi in enumerate(self.nr_orbitals):
            self.atom_orb_indices2[i,:noi]=self.atom_orb_indices[i]


    def get_number_of_transformations(self):
        '''
        Return the number of symmetry transformations in different directions.
        '''
        n = []
        r = self.atoms.get_symmetry_operation_ranges()
        for i in range(3):
            if r[i,0]==-nu.Inf:
                n.append(nu.Inf)
            else:
                n.append( int(round(r[i,1]-r[i,0]+1)) )
        return n
    
    def get_valences(self):
        """ Number of valence electrons for atoms. """
        return self.nr_of_valences

    def get_number_of_electrons(self):
        """ Total number of electrons in system. """
        return self.electrons

    def get_files(self):
        """ Return the dictionary for element input files. """
        return self.files

    def get_present(self):
        """ Return list of elements present. """
        return self.present

    def get_positions(self):
        return self.atoms.get_positions()/Bohr

    def get_center_of_mass(self):
        """ Return the center of mass. """
        return self.atoms.get_center_of_mass() / Bohr

    def get_atomic_numbers(self):
        """ Return the atomic numbers. """
        return self.atoms.get_atomic_numbers()

    def get_symbols(self):
        return self.symbols
    
    def get_N(self):
        """ Return the number of atoms. """
        return self.N

    def symbol(self,i):
        return self.symbols[i]

    def get_element(self,i):
        '''
        Return the element object of given atom.
        
        @param i: atom index or element symbol
        '''
        if isinstance(i,int):
            return self.elements[self.symbols[i]]
        else: 
            return self.elements[i]


    def orbitals(self,i=None,basis=False,atom=False,indices=False,number=False):
        """ Return data related to system orbitals.

        Parameters:
        -----------
        Default: return list of orbital-dictionaries for all atoms.
        i: return dict of orbitals for given atom  (dict incluldes 'atom','symbol','orbital',...)
        basis: return the data for basis orbital i.
        atom: return atom index for the basis orbital index i.
        indices: return orbital indices for atom i.
        number: return number of orbitals for atom i.
        """
        if i==None:
            return self.orb
        elif i!=None:
            if atom:
                return self.orbital_atoms[i]
            elif indices:
                return self.atom_orb_indices[i] #[orb['index'] for orb in self.atomorb[i]]
            elif number:
                return self.nr_orbitals[i] #len(self.atomorb[i])
            elif basis:
                return self.orb[i]
            else:
                return self.atomorb[i]


    def get_nr_orbitals(self):
        """ Total number of orbitals. """
        return self.norb
   
    
    def get_property_lists(self,lst=['i']):
        '''
        Return lists of atoms' given properties.
        
        @param lst: 'i'=index; 's'=symbol; 'no'=number of orbitals; 'o1'= first orbital 
        '''
        def get_list(p):
            if p=='i':      return range(self.N)
            elif p=='s':    return self.symbols
            elif p=='no':   return self.nr_orbitals
            elif p=='o1':   return self.first_orbitals
            else:
                raise NotImplementedError('Property not defined')
        l = [ get_list(item) for item in lst ]
        return zip(*l)
            
            
    def get_cube(self):
        """ Return the (orthorhombic) unit cell cube dimensions. """
        cell=abs(self.atoms.get_cell())
        e=1E-10
        if cell[0,1]>e or cell[0,2]>e or cell[1,2]>e:
            raise AssertionError('For cube the unit cell has to be orthorombic')
        return self.atoms.get_cell().diagonal()/Bohr
        

    def get_free_population(self, m):
        """ Return the population of the basis state m when the atom
        is isolated (a number between zero and two). """
        orb = self.orb[m]
        i_element = orb['atom']
        n_el = self.get_valences()[i_element]
        atomindex = orb['atomindex'] # atomindex states before this
        return max(0, min(2, n_el - 2*atomindex))
    
    
    def get_free_atoms_energy(self):
        '''
        Return the total of free atom energies for the system. 
        '''
        e = 0.0
        for s in self.symbols:
            e += self.elements[s].get_free_atom_energy()
        return e
    
    

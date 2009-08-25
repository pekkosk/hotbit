from element import Element
from ase import Atoms as ase_Atoms
from hotbit.atoms import Atoms
import numpy as nu
import hotbit.auxil as aux
import box.mix as mix
from numpy.linalg.linalg import norm
from ase.units import Hartree,Bohr
from os import environ
from weakref import proxy
from copy import copy, deepcopy

#from atoms import BravaisAtoms
#from atoms import WedgeAtoms




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
        #self.atoms = atoms.copy()
        self._update_atoms(atoms)
            #raise AssertionError('Given atoms object has to be box.Atoms, not ase.Atoms type.')
        #if isinstance(atoms,Atoms):
        #    self.atoms=atoms.copy()
        #else:
        #    self.atoms=Atoms(atoms)
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
        default_dir=environ.get('HOTBIT_PARAMETERS')
        self.files={}
        if elements==None or ('others' in elements and elements['others']=='default'):
            if elements!=None and 'others' in elements:
                elements.pop('others')
            for elem in self.present:
                self.files[elem]='%s/%s.elm' %(default_dir,elem)

        if elements!=None:
            self.files.update(elements)
        self._elements_initialization()
        self.solved={'ground state':None,'energy':None,'forces':None,'stress':None,'ebs':None,'ecoul':None}

    def __len__(self):
        return len(self.symbols)

    def __del__(self):
        pass
    

    def greetings(self):
        """ Return documentation for elements from .elm files. """
        txt='%i atoms, %i states, %.1f electrons (%.1f filled states)\n' %(self.N,self.norb,self.electrons,self.electrons/2)
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
        for i in range(3):
            assert r[i,0]<=r[i,1]
            if not self.atoms.get_pbc()[i]: assert r[i,0]==r[i,1]==0
            if r[i,0]==-nu.Inf:
                assert r[i,1]==nu.Inf
                r[i,:] = [-Mlarge,Mlarge]
            elif r[i,0]<-Mlarge:
                assert r[i,1]>=Mlarge
                r[i,:] = [-Mlarge,Mlarge]
            self.ranges.append( range(int(round(r[i,0])),int(round(r[i,1]))+1) )


    def calculation_required(self,atoms,quantities):
        """ Return True if quantities are solved for atoms.

        The quantities can be one or more of: 'ground state', 'energy', 'forces', and 'stress'.
        'ebs' or 'ecoul'
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
        self.Rn = [[self.nvector(r=i,ntuple=(0,0,0)) for i in xrange(self.N)]]
        self.Tn = [[self.nvector(r=i,ntuple=(0,0,0),lst='tensor') for i in xrange(self.N)]]
                
        cut2 = self.calc.ia.hscut**2
        self.ntuples = [(0,0,0)]
        # calculate the distances from unit cell 0 to ALL other possible; select meaningful
        for n1 in self.ranges[0]:
            for n2 in self.ranges[1]:
                for n3 in self.ranges[2]:
                    n = (n1,n2,n3)
                    if n==(0,0,0): continue
                    # check that any atom interacts with this unit cell
                    add = False    
                    R = [self.nvector(r=i,ntuple=n) for i in xrange(self.N)]
                    for i in xrange(self.N):
                        if add: break
                        for j in xrange(self.N):
                            dR = self.Rn[0][i]-R[j]
                            if dR[0]**2+dR[1]**2+dR[2]**2 <= cut2[i,j]: 
                                add = True
                                break
                    if add:
                        T = [self.nvector(r=i,ntuple=n,lst='tensor') for i in range(self.N)]
                        self.ntuples.append(n)
                        self.Rn.append(R)
                        self.Tn.append(T)
                        
        
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
        self.Rn = nu.array(self.Rn)
        self.Tn = nu.array(self.Tn)
        self.calc.stop_timing('geometry')
        
    def get_pbc(self):
        return self.atoms.get_pbc()

    def get_transforms(self):
        return self.ntuples

       
    def rotation_of_axes(self,n):
        '''
        Return the quantization axis rotation matrix for given symmetry operation.
         
        @param n: 3-tuple for transformation 
        '''
        return self.atoms.rotation_of_axes(n)
        
    
    def nvector(self,r,ntuple=(0,0,0),r0=[0,0,0],lst='vec'):
        '''
        Return position vector rn-r0, when r is operated by S(n)r=rn.
        
        @param r:   position (array) or atom index (integer) which is operated
        @param ntuple:   operate on r with S(n)
        @param r0:  if integer, use atom r0's position as r0
        @param l:   list of properties to return, 'vec'=vector, 'hat'=unit vector, 'norm'=norm
                    'tensor' = return the dyadic tensor (derivative) T(r,n)_ab = d rn_a/d r_b
        '''
        if not isinstance(lst,(list,tuple)):
            lst=[lst]
        assert not( r0!=[0,0,0] and 'tensor' in lst )
        if isinstance(r,int):  r=self.atoms.positions[r]
        if isinstance(r0,int): r0=self.atoms.positions[r0]
        vec=(self.atoms.transform(r,ntuple)-r0) / Bohr
        
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
            elif l=='tensor': 
                ret.append( self.atoms.tensor(r,ntuple) )
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
        self.norb=0
        for i,symb in enumerate(self.symbols):
            el=self.elements[symb]
            atomorb=[]
            self.first_orbitals.append(self.norb)
            for k,(ao,e) in enumerate(zip(el.get_orbital_types(),el.get_onsite_energies())):
                self.norb+=1
                Rnl=el.get_Rnl_function(ao)
                atomorb.append({'atom':i,'symbol':symb,'orbital':ao,\
                                'index':self.norb-1,'energy':e,'atomindex':k,'Rnl':Rnl})
            self.nr_orbitals.append(len(atomorb))
            self.atomorb.append(atomorb)
            self.orb.extend(atomorb)

        # number of valence electrons on atoms and the total number of electrons
        self.nr_of_valences=nu.array( [self.elements[symb].get_valence_number() for symb in self.symbols] )
        self.electrons=self.get_valences().sum()-self.charge

        # set list of element pairs
        self.el_pair_list=[]
        for i,si in enumerate(self.present):
            for j,sj in enumerate(self.present):
                self.el_pair_list.append([si,sj])

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

    def get_element_pairs(self):
        """ Return list of element pairs. """
        return self.el_pair_list

    def orbitals(self,i=None,indices=False,number=False):
        """ Return data related to system orbitals.

        Parameters:
        -----------
        Default: return list of orbital-dictionaries for all atoms.
        i: return dict of orbitals for given atom  (dict incluldes 'atom','symbol','orbital',...)
        indices: return orbital indices for atom i.
        number: return number of orbitals for atom i.
        """
        if i==None:
            return self.orb
        elif i!=None:
            if indices:
                return self.atom_orb_indices[i] #[orb['index'] for orb in self.atomorb[i]]
            elif number:
                return self.nr_orbitals[i] #len(self.atomorb[i])
            else:
                return self.atomorb[i]

    def get_nr_orbitals(self):
        """ Total number of orbitals. """
        return self.norb

    def get_ia_atom_pairs(self,get):
        # TODO: get rid of this!
        """
        List interacting atom pairs with desired properties:

        get : list of strings to get.
            Options: 'i','j','si','sj','r','rhat','dist','ri','rj','o1i','o1j','noi','noj','oi','oj'
            ( 'i'=atom index, 'si'=symbol, 'r'=vector i-->j, 'rhat'=unit vector, 'dist'=distance
            'ri'=position of i, 'o1i'=index of first orbital, 'noi'=number of orbitals,
            'oi'=orbital indices for i, and some anagolously for j.)
        """
        lst=[]
        for item in get:
            lst.append(self.ia_pairs[item])
        return zip(*lst)
    
    
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
            

    def get_atom_pairs(self,mode=1,cut=None,include=None):
        # TODO get rid of this method !!
        """
        List atom pairs with distance cutoff / only given element pairs.
        Different modes:
        1: return (i,j)
        2: return (i,si,j,sj)
        3: return (i,si,j,sj,dist)
        4: return (i,si,j,sj,dist,r)             (r from i to j)
        5: return (i,si,j,sj,dist,r,rhat)        (    -"-      )
        6: return (i,si,ri,j,sj,rj,dist,r,rhat)  (    -"-      )
        7: return (i,o1i,noi,j,o1j,noj) (o1i=first orbital index, noi=number of orbitals)
        """
        if include==None:
            included=self.get_element_pairs()
        else:
            included=include

        lst=[]
        if cut!=None:
            assert cut<=self.cut
        for (i,j) in self.ia_pairs:
                si, sj=self.symbols[i], self.symbols[j]
                if not ([si,sj] in included or [sj,si] in included):
                    continue
                d=self.geometry[(i,j)]['dist']
                if cut==None or (cut!=None and d<cut):
                    if mode==4 or mode==5 or mode==6:
                        rvec=self.geometry[(i,j)]['r']
                    if mode==1:
                        lst.append((i,j))
                    elif mode==2:
                        lst.append((i,si,j,sj))
                    elif mode==3:
                        lst.append((i,si,j,sj,d))
                    elif mode==4:
                        lst.append((i,si,j,sj,d,rvec))
                    elif mode==5:
                        lst.append((i,si,j,sj,d,rvec,rvec/norm(rvec)))
                    elif mode==6:
                        ri=self.vector(i)
                        rj=self.vector(j)
                        lst.append((i,si,ri,j,sj,rj,d,rvec,rvec/norm(rvec)))
                    elif mode==7:
                        noi, noj=self.nr_orbitals[i], self.nr_orbitals[j]
                        o1i, o1j=self.first_orbitals[i], self.first_orbitals[j]
                        noi, noj=self.nr_orbitals[i], self.nr_orbitals[j]
                        lst.append((i,o1i,noi,j,o1j,noj))
        return lst


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
from element import Element
from ase import Atoms
import numpy as nu
import hotbit.auxil as aux
import box.mix as mix
from numpy.linalg.linalg import norm
from ase.units import Hartree,Bohr
from os import environ
from weakref import proxy

from atoms import BravaisAtoms




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

        if not isinstance(atoms, Atoms):
            raise AssertionError('Given atoms object has to be ase.Atoms type.')
        self.atoms = atoms.copy()
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
        self._initialize()
        self.solved={'ground state':None,'energy':None,'forces':None,'stress':None}


    def init2(self):
        # TODO: these have to be divided into subroutines
        self.nranges=[(-1,1),(-1,1),(-1,1)]
#        self.nranges=[(0,0),(0,0),(0,0)]
        self.nlists= [range(n[0],n[1]+1) for n in self.nranges]
        self.M=[len(l) for l in self.nlists]
#        self.ntuples=[] # n-tuples
#        for n1 in self.nlists[0]:
#            for n2 in self.nlists[1]:
#                for n3 in self.nlists[2]:
#                    self.ntuples.append((n1,n2,n3))        

        self.ntuples=[] # n-tuples
        rang=[]
        for i in range(3):
            if self.atoms.pbc[0]:
                rang.append([0,1,-1])
            else:
                rang.append([0])
            
        for n1 in rang[0]:
            for n2 in rang[1]:
                for n3 in rang[2]:
                    self.ntuples.append((n1,n2,n3))
                    
#        self.copies = len(self.ntuples)
        self.Rn = nu.zeros((self.N,len(self.ntuples),3))
        self.Tn = nu.zeros((self.N,len(self.ntuples),3,3))
        for i in range(self.N):
            for n,nt in enumerate(self.ntuples):
                self.Rn[i,n,:] = self.nvector(r=i,ntuple=nt)
                self.Tn[i,n,:,:] = self.nvector(r=i,ntuple=nt,lst='dtensor')
            

    def get_transforms(self):
        return self.ntuples
       
    
    def nvector(self,r,ntuple=(0,0,0),r0=[0,0,0],lst='vec'):
        '''
        Return position vector rn-r0, when r is operated by S(n)r=rn.
        
        @param r:   position (array) or atom index (integer) which is operated
        @param ntuple:   operate on r with S(n)
        @param r0:  if integer, use atom r0's position as r0
        @param l:   list of properties to return, 'vec'=vector, 'hat'=unit vector, 'norm'=norm
                    'dtensor' = return the dyadic tensor (derivative) T(r,n)_ab = d rn_a/d r_b
        '''
        self.calc.start_timing('nvec')
        if not isinstance(lst,(list,tuple)):
            lst=[lst]
        assert not( r0!=[0,0,0] and 'dtensor' in lst )
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
            elif l=='dtensor': 
                ret.append( self.atoms.dtensor(r,ntuple) )
            else:
                raise AssertionError('Keyword %s not defined' %l)
        
        self.calc.stop_timing('nvec')
        if len(ret)==1: 
            return ret[0]
        else: 
            return ret

    def __del__(self):
        pass

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


    def set_name(self, name):
        self.name = name


    def __len__(self):
        return len(self.symbols)


    def _initialize(self):
        """ Initialize element objects, orbital tables etc. Do it only once."""
        self.elements={}
        for symb in self.present:
            if symb not in self.files:
                raise KeyError('Element file for %s was not defined.' %symb)
            if type(self.files[symb])==type(''):
                self.elements[symb]=Element(self.files[symb])
            else:
                self.elements[symb]=self.files[symb]

        self.orb=[]
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



    def update_geometry(self):
        # TODO: this method will become useless
        """ Update all properties related to geometry (calculate once/geometry)

        (analogously for j): i=atom index, si=symbol, r=vector i-->j, rhat=unit vector,
        ri=position of i, o1i=index of first orbital, noi=number of orbitals, oi=orbital indices of i
        """
        self.calc.start_timing('geometry')
        self.geometry={}
        properties=['i','j','si','sj','r','dist','rhat','ri','rj','o1i','o1j','noi','noj','oi','oj']
        self.ia_pairs={}
        for property in properties:
            self.ia_pairs[property]=[]
        self.cut=self.calc.ia.get_cutoffs()

        for i in range(self.N):
            for j in range(i,self.N):
                dist=self.distance(i,j)
                r=self.vector(i,j)
                rhat=r/dist
                self.geometry[(i,j)]={'dist':dist,'r':r,'rhat':rhat}
                if i!=j:
                    self.geometry[(j,i)]={'dist':dist,'r':-r,'rhat':-rhat}
                si, sj=self.symbols[i], self.symbols[j]
                if dist<=self.cut[si+sj]:
                    ri, rj=self.vector(i), self.vector(j)
                    o1i, o1j=self.first_orbitals[i], self.first_orbitals[j]
                    noi, noj=self.nr_orbitals[i], self.nr_orbitals[j]
                    oi, oj=self.atom_orb_indices[i], self.atom_orb_indices[j]
                    for key,value in zip(properties,[i,j,si,sj,r,dist,rhat,ri,rj,o1i,o1j,noi,noj,oi,oj]):
                        self.ia_pairs[key].append(value)

        # setup a table which orbitals are interacting
        self.ia_orbitals=nu.zeros((self.norb,self.norb),int)
        self.nr_ia_orbitals=nu.zeros((self.norb,),int)
        for oi, oj, noj in self.get_ia_atom_pairs(['oi','oj','noj']):
            # all orbitals oi and oj interact with each other
            for orbi in oi:
                ni=self.nr_ia_orbitals[orbi]
                self.ia_orbitals[orbi,ni:ni+noj]=oj
                self.nr_ia_orbitals[orbi]+=noj
        self.calc.stop_timing('geometry')


    def calculation_required(self,atoms,quantities):
        """ Return True if quantities are solved for atoms.

        The quantities can be one or more of: 'ground state', 'energy', 'forces', and 'stress'.
        """
        if not isinstance(quantities,(list,tuple)):
            quantities=[quantities]

        if type(self.solved['ground state']) == type(None):
            return True

        # check that all quantities have been solved for identical atoms
        for quantity in quantities:
            solved_atoms=self.solved[quantity]
            if type(solved_atoms) == type(None):
                return True
            elif not atoms == solved_atoms:
                return True
        return False


    def set_atoms(self,atoms):
        """ Set the atoms object ready for calculations. """
        if type(atoms) == type(None) or not atoms == self.atoms:
#            self.atoms=Atoms(atoms)
            self.atoms=BravaisAtoms(atoms)
            self.update_geometry()
            self.init2()


    def set_solved(self,quantities):
        """ Set quantities solved for current atoms. """
        if not isinstance(quantities,(list,tuple)):
            quantities=[quantities]
        for quantity in quantities:
            # FIXME: more careful bookkeeping for solved stuff 
            self.solved[quantity]=BravaisAtoms(self.atoms)


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


    def vector(self, ri, rj=None, mic=True):
        """
        Return the vector ri or rj-ri, normally within the minimum image convention (mic).
        
        parameters:
        -----------
        ri: initial point (atom index or position)
        rj: final point (atom index or position)
        """
        if not mic:
            raise NotImplementedError('Why not use mic?')
        R = self.atoms.get_positions() / Bohr
        if type(ri) == type(1):
            Ri = R[ri].copy()
        else:
            Ri = ri.copy()
        if rj == None:
            Rj = None
        elif type(rj) == type(1):
            Rj = R[rj].copy()
        else:
            Rj = rj

        if Rj == None:
            rij = Ri
        else:
            rij = Rj - Ri

        L = self.get_cell_axes()
        for a in range(3):
            if self.atoms.pbc[a]:
                rij[a] = nu.mod(rij[a]+0.5*L[a],L[a]) - 0.5*L[a]
        return nu.array(rij)


    def distance(self,ri,rj=None,mic=True):
        """
        Return the length of |ri| or |rj-ri|, normally within the minimum image convention.
        """
        if rj == None:
            return mix.norm(self.vector(ri,mic))
        else:
            return mix.norm(self.vector(ri,rj,mic))


    def distance_of_elements(self, si, sj, mode='minimum'):
        """ Returns the closest distance between two elements. """
        if type(si) != str or type(sj) != str:
            raise AssertionError('Wrong types of parameters')
        found_pair = False
        distances = []
        for i, a in enumerate(self.atoms):
            for j, b in enumerate(self.atoms):
                if i == j:
                    pass
                else:
                    if (a.symbol == si and b.symbol == sj):
                        found_pair = True
                        d = self.distance(i, j)
                        distances.append(d)
        if found_pair:
            if mode == 'minimum':
                return min(distances)
            else:
                raise NotImplementedError('Only mode=minimum works')
        else:
            return None


    def get_cell_axes(self):
        """
        Lengths of the unit cell axes (currently for orthorhombic cell).
        """
        x,y,z = self.atoms.get_cell()
        assert (nu.dot(x,y) == 0 and nu.dot(x,z) == 0 and nu.dot(y,z) == 0)
        return nu.diag(self.atoms.get_cell()) / Bohr


    def get_box_lengths(self):
        """ Return lengths of box sizes in orthorombic box. """
        c=self.atoms.get_cell()/Bohr
        return c[0,0],c[1,1],c[2,2]


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

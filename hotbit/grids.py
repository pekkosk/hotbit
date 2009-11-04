import numpy as nu
from weakref import proxy
from ase.units import Bohr
from hotbit.analysis.wavefunctions import angular
from math import sqrt
from scipy.special import erf


class Grids:
    def __init__(self,calc):
        self.calc = proxy(calc)
        self.spacing = None
        self.grid = None
        self.ng = None
        self.gbasis = None
        self.el=proxy(self.calc.el)
        self.Rcut = 5.0/Bohr  # cutoff for atomic radial parts
        self.Rcore = 0.25/Bohr # the core-region

        
    def cut_core(self,d):
        """ Cut the core wiggles with this screening function. 
        
        Wiggles have to be cut at least within the length scale of the grid spacing.
        """
        scale = self.spacing*5 
        return erf(d/scale)
        
        
    def make_grid(self,spacing=0.2/Bohr,pad=False):
        """
        Make grid.
        
        @param spacing: grid spacing in Bohrs
        """
        self.calc.start_timing('make grid')
        if any(self.calc.el.atoms.get_pbc()):
            raise AssertionError('Grid stuff not implemented for periodic systems yet.')
        if self.spacing==spacing:
            return
        else:
            self.grid, self.ng, self.gbasis = None, None, None
        self.norb = self.calc.get_number_of_bands()
        self.spacing=spacing
        self.L = self.el.get_cube()        
        self.N = nu.ceil( self.L/spacing )
        self.dr = self.L/self.N
        if not pad:
            self.L -= self.dr # this is due to the cube-format
        self.gd1=[]
        for i in range(3):
            self.gd1.append( nu.linspace(0,self.L[i],self.N[i]) )
        
        self.dV = nu.prod(self.dr)
        self.ng = nu.prod(self.N)
        self.grid = []
        for i,x in enumerate(self.gd1[0]):
            for j,y in enumerate(self.gd1[1]):
                for k,z in enumerate(self.gd1[2]):
                    self.grid.append(nu.array([x,y,z]))
        self.calc.stop_timing('make grid') 
        
    
    def basis_orbital_to_grid(self,b):
        """ 
        Put (real) basis orbital into grid. 
        
        @param b: index of basis orbital.
        
        Normalize the orbital by scaling wf within the core region.
        """
        self.calc.start_timing('to grid')
        if self.gbasis!=None:
            return self.gbasis[b]
        
        basis = nu.zeros(self.ng)
        atom = self.el.orbitals(b,atom=True)
        orb = self.el.orbitals(b,basis=True)
        nvector = self.el.nvector
        type, Rnl = orb['orbital'], orb['Rnl']
        core = []
        for ig,r in enumerate(self.grid):    
            dr = nvector(r,r0=atom)
            d = sqrt(dr[0]**2+dr[1]**2+dr[2]**2)
            if d>self.Rcut: 
                continue
            else:
                basis[ig] = Rnl(d)*angular(dr,type)*self.cut_core(d)
                if d<self.Rcore: core.append(ig)
        assert len(core)>1
        # normalize basis function; correct wf WITHIN CORE region 
        # to get normalization to one. 
        nall = sum(basis**2)*self.dV
        ncore = sum(basis[core]**2)*self.dV
        basis[core] = basis[core] * sqrt((1-nall+ncore)/(ncore))  
        basis.shape = tuple(self.N)
        self.calc.stop_timing('to grid')
        return basis
    
    
    def whole_basis_to_grid(self):
        """ Put all basis functions to grid. """
        self.calc.start_timing('basis2grid')
        gbasis = [ self.basis_orbital_to_grid(b) for b in range(self.norb) ]
        self.gbasis = nu.array(gbasis)
        self.calc.stop_timing('basis2grid')
            
    
    def get_wf(self,i,k=0):
        """ 
        Put (complex) wave function into grid. 
        
        @param i: band index
        @param k: k-point index
        The wave function is NOT normalized.
        """
        if self.gbasis == None:
            self.whole_basis_to_grid()
        wf = self.calc.st.wf
        # TODO: this can be done more efficiently using numpy
        gwf = nu.zeros_like( self.gbasis[0] ) #TODO: complex
        for b in range(self.norb):
            gwf = gwf + self.gbasis[b]*wf[k,b,i] #TODO fix this 
        return gwf       
    
    
    def get_wf_density(self,i,k=0):
        """
        Return normalized density from given eigenfunction.
        
        @param i: band index
        @param k: k-point index
        
        Density is normalized by rescaling. 
        """
        wf = self.get_wf(i,k)
        dens = wf*wf.conjugate()
        return dens/(dens.sum()*self.dV)
    
    
    
    def get_density(self):
        """ 
        Return valence electron density on grid.
        
        Integrates exactly to the total number of valence electrons. 
        """
        k=0        
        dens = nu.zeros(tuple(self.N))
        for i,f in enumerate(self.calc.st.f[k,:]):
            if f<1E-6: break
            dens = dens + self.get_wf_density(i,k)*float(f) #TODO: k-point weights 
        assert all(abs(dens.flatten().imag)<1E-10)
        assert abs(dens.flatten().sum()*self.dV-self.calc.el.get_number_of_electrons())<1E-9
        return dens.real
        
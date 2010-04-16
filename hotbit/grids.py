import numpy as nu
from weakref import proxy
from ase.units import Bohr
from hotbit.analysis.wavefunctions import angular
from math import sqrt
from scipy.special import erf


class Grids:
    def __init__(self,calc,spacing,pad):
        """
        Initialize grids.
        
        parameters:
        ===========
        spacing: grid spacing in Angstroms
        pad:     True for padded edges (grid points at opposite edges 
                 have the same value)
        """
        self.calc = proxy(calc)
        self.el=proxy(self.calc.el)
        self.spacing = spacing/Bohr
#        self.grid = None
#        self.ng = None
#        self.gbasis = None
  
#        self.Rcut = 5.0/Bohr  # cutoff for atomic radial parts
 #       self.Rcore = 0.3/Bohr # the core-region
#        self.maxspacing = self.Rcore/sqrt(3*0.5**2)
 #       self._make_grid(spacing,pad)
             
        self.calc.start_timing('init grid')
  #      if nu.any(self.calc.el.atoms.get_pbc()):
  #          raise AssertionError('Grid stuff not implemented for periodic systems yet.')
   #     if spacing>self.maxspacing:
   #         raise AssertionError('Grid spacing must be smaller than %.3f Angstroms' %(self.maxspacing*Bohr))
   #     if self.spacing==spacing:
   #         return
   #     else:
   #         self.grid, self.ng, self.gbasis = None, None, None
    #    self.norb = self.calc.get_number_of_bands()
   #     self.spacing=spacing
        self.L = self.el.get_cube()        
        self.N = nu.array( nu.ceil(self.L/spacing),int )
        self.dr = self.L/self.N
        if not pad:
            self.L -= self.dr # this is due to the .cube-format
        self.grid=[]
        for i in range(3):
            self.grid.append( nu.linspace(0,self.L[i],self.N[i]) )
        
        
        self.dV = nu.prod(self.dr)
        self.ng = nu.prod(self.N)
        #self.grid = []
        #for i,x in enumerate(self.gd1[0]):
        #    for j,y in enumerate(self.gd1[1]):
        #        for k,z in enumerate(self.gd1[2]):
        #            self.grid.append(nu.array([x,y,z]))
        
        #self.grid = nu.zeros((self.N[0],self.N[1],self.N[2],3))
        #for i,x in enumerate(self.gd1[0]):
        #    for j,y in enumerate(self.gd1[1]):
        #        for k,z in enumerate(self.gd1[2]):
        #            self.grid[i,j,k] = (x,y,z)        
                    
        # partial grids
        self.cutoff = 5.0/Bohr # cutoff for atomi wfs
        self.pN = []
        self.pL = []
        for i in range(3):
            N = int( nu.ceil(self.cutoff*2/self.dr[i]) )
            if nu.mod(N,2)==0: 
                N+=1
            self.pN.append(N)
            self.pL.append(N*self.dr[i])
        
        self.pL = nu.array(self.pL)
        self.pgrid = []
        for i in range(3):
            #self.pgrid.append( nu.linspace(0,self.L[i],self.N[i]) )
            self.pgrid.append( [p*self.dr[i]-self.pL[i]/2 for p in xrange(self.pN[i])] )
            
        #self.pgrid = nu.zeros((self.pN[0],self.pN[1],self.pN[2],3))
        #for i in range(self.pN[0]):
        #    for j in range(self.pN[1]):
        #        for k in range(self.pN[2]):
        #            self.pgrid[i,j,k] = (i,j,k)*self.dr - self.pL/2  
                    
        self._all_basis_orbitals_to_partial_grid()             
        self.calc.stop_timing('init grid') 
        
    
    def _basis_orbital_to_partial_grid(self,symbol,otype):
        """
        Return atomic orbital on a grid.
        
        parameters:
        ===========
        symbol:    atom symbol
        otype:     orbital type ('s','px','py',...)
        """
        el = self.calc.el.elements[symbol]
        Rnl = el.get_Rnl_function(otype) 
        range = el.get_wf_range(otype,fractional_limit=1E-5)
        
        wf = nu.zeros(self.pN)
        for i in xrange(self.pN[0]):
            for j in xrange(self.pN[1]):
                for k in xrange(self.pN[2]):
                    r = nu.array( [self.pgrid[0][i],self.pgrid[1][j],self.pgrid[2][k]] )
                    d = sqrt( r[0]**2+r[1]**2+r[2]**2 )
                    if d>self.cutoff or d>range:
                        continue
                    else:
                        rnl = Rnl(d)
                        if abs(rnl)>1E-7:
                            wf[i,j,k] = rnl*angular(r,otype)
        return wf     
    
    
    def _all_basis_orbitals_to_partial_grid(self):
        """
        Put all basis orbital types into grid 
        """
        self.calc.start_timing('orbitals to grid') 
        self.basis={}
        for symb in self.el.present:
            el = self.el.elements[symb]
            symbol = el.get_symbol()
            self.basis[symbol] = {}
            for otype in el.get_orbital_types():
                wf = self._basis_orbital_to_partial_grid(symbol,otype)
                self.basis[symbol][otype] = wf
        self.calc.stop_timing('orbitals to grid')
    

    def _basis_orbital_to_grid(self,I,otype):
        """
        Put basis orbital into the full grid.
        """
        symbol = self.el.symbols[I]
        
        ri = self.el.nvector(I)
        #position of atom I => grid point N
        N = nu.array( nu.round(ri/self.dr),int ) 
        dn = []
        a,b = nu.zeros((3,2),int),nu.zeros((3,2),int)
        for i in range(3):
            mn = N[i]-(self.pN[i]-1)/2        
            mx = N[i]+(self.pN[i]-1)/2
            # first is index, second is the NUMBER of items (not index)
            a[i] =( max(0,mn),min(self.N[i],mx+1) )
            b[i] =( max(0,-mn),min(self.pN[i],self.pN[i]+self.N[i]-(mx+1)) )
            #print i,mn,mx,a[i],b[i],a[i,1]-a[i,0], b[i,1]-b[i,0]
        wf = nu.zeros(self.N)
        pwf = self.basis[symbol][otype]
#        print '--'*10
#        print self.N
#        print self.pN
#        print N
#        print a[0,0],a[0,1],a[1,0],a[1,1],a[2,0],a[2,1]
#        print b[0,0],b[0,1],b[1,0],b[1,1],b[2,0],b[2,1]
        wf[a[0,0]:a[0,1],a[1,0]:a[1,1],a[2,0]:a[2,1]] = pwf[b[0,0]:b[0,1],b[1,0]:b[1,1],b[2,0]:b[2,1]]           
        return wf
        
    
    def get_wf(self,a,k=0):
        """ 
        ff
        """        
        wf = self.calc.st.wf
        # TODO: this can be done more efficiently using numpy
        gwf = nu.zeros(self.N) #TODO: complex
        
        for mu,orb in enumerate(self.el.orbitals()):
            I,otype = orb['atom'],orb['orbital']
            gwf += wf[k,a,mu]*self._basis_orbital_to_grid(I,otype)
        return gwf
                 
            
        
        
    def basis_orbital_to_grid(self,b):
        """ 
        Put (real) basis orbital into grid. 
        
        @param b: index of basis orbital.
        
        Normalize the orbital by scaling wf within the core region.
        """
        self.calc.start_timing('to grid')
        if self.gbasis!=None:
            return self.gbasis[b]
        #aprint 'basis',b
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
                basis[ig] = Rnl(d)*angular(dr,type)
                if d<self.Rcore: core.append(ig)
        # normalize basis function; correct wf WITHIN CORE region 
        # to get normalization to one. 
        nall = sum(basis**2)*self.dV
        if nall<0.1:
            raise AssertionError('Some element probably does not have radial functions.')
        ncore = sum(basis[core]**2)*self.dV
        if ncore<1E-10:
            raise AssertionError('Wf core correction not possible; spacing is probably too large.')
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
            
    
    def get_wf_old(self,i,k=0):
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
        assert nu.all(abs(dens.flatten().imag)<1E-10)
        assert abs(dens.flatten().sum()*self.dV-self.calc.el.get_number_of_electrons())<1E-9
        return dens.real
        
import numpy as np
from weakref import proxy
from ase.units import Bohr
from hotbit.analysis.wavefunctions import angular
from math import sqrt
from scipy.special import erf


class Grids:
    def __init__(self,calc,h,cutoff=3.0):
        """
        Initialize grids.
        
        parameters:
        ===========
        h:       grid spacing in Angstroms
        pad:     True for padded edges (grid points at opposite edges 
                 have the same value)
        cutoff:  cutoff for atomic orbitals
        """
        self.calc = proxy(calc)
        self.el = proxy(self.calc.el)
        self.h = h/Bohr
        self.cutoff = cutoff/Bohr # cutoff for atomic orbitals
        self.clip = 100.0 #clip wfs             
        self.calc.start_timing('init grid')
        
        
        self.L = self.el.get_cube()        
        self.N = np.array( np.round(self.L/h),int )
        self.dr = self.L/(self.N-1)

        self.grid=[]
        for i in range(3):
            self.grid.append( np.linspace(0,self.L[i],self.N[i]) )

        
        self.dV = np.prod(self.dr)
        self.ng = np.prod(self.N)

        # partial grids (where atomic orbitals are first put)        
        self.pN = []
        self.pL = []
        for i in range(3):
            N = int( np.round(self.cutoff*2/self.dr[i]) )
            if np.mod(N,2)==0: 
                N+=1
            self.pN.append(N)
            self.pL.append((N-1)*self.dr[i])
        
        self.pL = np.array(self.pL)
        self.pgrid = []
        for i in range(3):
            self.pgrid.append( [p*self.dr[i]-self.pL[i]/2 for p in xrange(self.pN[i])] )

                    
        self._all_basis_orbitals_to_partial_grid()             
        self.calc.stop_timing('init grid') 
        
        
    def _return_array(self,a,pad):
        a=a.copy()
        a = a.real.clip(-self.clip,self.clip) + 1j*a.imag.clip(-self.clip,self.clip)
        if np.all(abs(a.imag)<1E-10):
            a=a.real
        if pad:
            return a
        else:
            return a[:-1,:-1,:-1]
        
    
    def _atomic_orbital_to_partial_grid(self,symbol,otype):
        """
        Put atomic orbital on a grid.
        
        It looks like this; atom is right in the middle of the grid
        (number of grid points is odd in each direction)
        
        |----------v----------|
        |     |    |    |     |
        |---------------------|
        |     |    |    |     |
        |----------X----------|
        |     |    |    |     |
        |---------------------|
        |     |    |    |     |
        |----------^-----------
        
        
        parameters:
        ===========
        symbol:    atom symbol
        otype:     orbital type ('s','px','py',...)
        """
        el = self.calc.el.elements[symbol]
        Rnl = el.get_Rnl_function(otype) 
        range = el.get_wf_range(otype,fractional_limit=1E-5)
        
        wf = np.zeros(self.pN)
        for i in xrange(self.pN[0]):
            for j in xrange(self.pN[1]):
                for k in xrange(self.pN[2]):
                    r = np.array( [self.pgrid[0][i],self.pgrid[1][j],self.pgrid[2][k]] )
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
        Put all atomic orbitals into partial grids -> e.g. patomic['C']['px']
        """
        self.calc.start_timing('atomic orbitals to grid') 
        self.patomic={}
        for symb in self.el.present:
            el = self.el.elements[symb]
            symbol = el.get_symbol()
            self.patomic[symbol] = {}
            for otype in el.get_orbital_types():
                wf = self._atomic_orbital_to_partial_grid(symbol,otype)
                self.patomic[symbol][otype] = wf
        self.calc.stop_timing('atomic orbitals to grid')
        
        
    
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
        symbol = self.el.symbols[I]
        pwf = self.patomic[symbol][otype]
        wf = np.zeros(self.N)
        phases = self.calc.ia.get_phases()[:,k]
        for ni,n in enumerate(self.el.ntuples):
            ri = self.el.nvector(I,n)  
            inside = True
            # check first roughly that basis reaches the inside of cell
            for i in range(3):
                if ri[i]<-self.cutoff or self.L[i]+self.cutoff<ri[i]:
                    inside = False
            if not inside:
                continue
            #position of atom I => grid point N
            N = np.array( np.round(ri/self.dr),int )
            a,b = [],[]
            
            for i in range(3):
                lo = N[i]-(self.pN[i]-1)/2        
                hi = N[i]+(self.pN[i]-1)/2
                # these are indices
                a1,a2 = max(0,lo), min(self.N[i]-1,hi)
                an = a2-a1
                b1 = max(0,-lo)
                b2 = b1 + an
                a.append( slice(a1,a2+1) )
                b.append( slice(b1,b2+1) )
                if b1>b2 or a1>a2:
                    inside = False
            if inside:
                wf[a[0],a[1],a[2]] = wf[a[0],a[1],a[2]] + phases[ni]*pwf[b[0],b[1],b[2]]
        return self._return_array(wf,pad) 

    
    
    def get_grid_wf(self,a,k=0,pad=True):
        """ 
        Return eigenfunction on a grid.
        
        parameters:
        ===========
        a:     state (band) index
        k:     k-vector index
        pad:   padded edges 
        """        
        wf = self.calc.st.wf
        
        if np.all( abs(self.calc.ia.get_phases()[:,k]-1)<1E-10 ):
            gwf = np.zeros(self.N) #TODO: complex
        else: 
            gwf = np.zeros(self.N,complex)
        
        for mu,orb in enumerate(self.el.orbitals()):
            I,otype = orb['atom'],orb['orbital']
            if abs(wf[k,a,mu])**2>1E-13:
                gwf += wf[k,a,mu]*self.get_grid_basis_orbital(I,otype,k)
        return self._return_array(gwf,pad)
            
    
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
        wf = self.get_grid_wf(a,k,pad=True)
        return self._return_array(wf*wf.conjugate(),pad) 
    
    
    def get_grid_density(self,pad):
        """ 
        Return electron density on grid.
        
        Do not perform accurate analysis on this density.
        Integrated density differs from the total number of electrons.
        Bader analysis will be inaccurate.
        
        parameters:
        pad:      padded edges
        """
        rho = np.zeros(tuple(self.N))
        for k,wk in enumerate(self.calc.st.k):
            for a,f in enumerate(self.calc.st.f[k,:]):
                if f<1E-6: break
                rho = rho + wk*self.get_grid_wf_density(a,k,pad=True)
        return self._return_array(rho,pad)    
        
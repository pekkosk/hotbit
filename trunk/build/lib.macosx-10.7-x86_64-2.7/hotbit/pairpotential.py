# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

import numpy as np
from ase.units import Hartree, Bohr
from weakref import proxy

def ij2s(i,j):
    return '%04i%04i' %(i,j)
    
    

class PairPotential:
    def __init__(self,calc):
        self.calc=proxy(calc)
        self.N=calc.el.get_N()
        self.ex = False
        self.comment = ''            
        self.atomv = {}
        self.elemv = {}
        self.rcut=0.0
        self.symbols = [s[0] for s in self.calc.el.get_property_lists(['s'])]
        self.greet = False
        
        
    def exists(self):
        """ Return whether any potentials exist. """
        return self.ex
        

    def greetings(self):
        """ Return the repulsion documentations. """
        print>> self.calc.txt,'Pair potentials:\n'+self.comment
    

    def add_pair_potential(self,i,j,v,eVA=True):
        """
        Add pair interaction potential function for elements or atoms.
        
        parameters:
        ===========
        i,j:    * atom indices, if integers (0,1,2,...)
                * elements, if strings ('C','H',...)
        v:      Pair potential function. 
                Only one potential per element and atom pair allowed. 
                Syntax:  v(r,der=0), v(r=None) returning the
                interaction range in Bohr or Angstrom.
        eVA:    True for v using  eV and Angstrom
                False for v using Hartree and Bohr
        """
        # construct a function that works in atomic units
        if eVA:
            def v2(r,der=0):
                if r==None:
                    return v(None)/Bohr
                if der==0:
                    return v(r/Bohr,der)/Hartree
                elif der==1:
                    return v(r/Bohr,der)*Bohr/Hartree
        else:
            v2 = v           
            
        if isinstance(i,int):
            self.comment += '  Atom pair %i-%i\n' %(i,j)
            self.atomv[ij2s(i,j)] = v2
            self.atomv[ij2s(j,i)] = v2
        else:
            self.comment += '  Element pair %s-%s\n' %(i,j)
            self.elemv[i+j]=v2
            self.elemv[j+i]=v2
        self.rcut = max(self.rcut, v2(None))
        self.ex = True
        
        
    def _get_full_v(self,i,j):
        """ Return the full interaction potential between atoms i and j."""
        ij = ij2s(i,j)
        sij = self.symbols[i]+self.symbols[j]
        
        if ij in self.atomv and sij in self.elemv:
            def v(r,der=0): 
                return self.atomv[ij](r,der) + self.elemv[sij](r,der)
        elif ij in self.atomv:
            def v(r,der=0): 
                return self.atomv[ij](r,der)
        elif sij in self.elemv:
            def v(r,der=0): 
                return self.elemv[sij](r,der)
        else:
            def v(r,der=0):
                return 0.0
        return v      


    def get_energy(self):
        """ Return the energy in pair potentials (in eV). """
        if not self.ex:
            return 0.0
        self.calc.start_timing('e_pp')
        if not self.greet:
            self.greetings()
            self.greet = True            
        
        lst=self.calc.el.get_property_lists(['i','s'])
        Rijn = self.calc.el.rijn
        epp=0.0
        for i,si in lst: 
            for j,sj in lst[i:]:
                for n,rijn in enumerate(Rijn[:,i,j]): 
                    if i==j and n==0: continue
                    d = np.sqrt( rijn[0]**2+rijn[1]**2+rijn[2]**2 )
                    if d>self.rcut: continue
                    v = self._get_full_v(i,j)                 
                    if i==j:
                        epp+=0.5*v(d)
                    else:
                        epp+=v(d)
        self.calc.stop_timing('e_pp') 
        return epp*Hartree
    
    
    def get_pair_energy(self,i,j):
        """
        Return the pair repulsion energy of given atom pair (in Hartree)
        
        parameters:
        ===========
        i,j:     atom indices
        """
        Rijn = self.calc.el.rijn
        epp=0.0
        for n,rijn in enumerate(Rijn[:,i,j]): 
            if i==j and n==0: continue
            d = np.sqrt( rijn[0]**2+rijn[1]**2+rijn[2]**2 )
            if d>self.rcut: continue
            v = self._get_full_v(i,j)
            if i==j:
                epp += 0.5*v(d)
            else:
                epp += v(d)
        return epp
    

    def get_forces(self):
        """ 
        Return pair potential forces in atomic units.
        
        F_i = sum_(j,n) V'(rijn) rijn/dijn, with rijn = r_j^n -r_i and dijn=|rijn|
        """
        f=np.zeros((self.N,3))
        if not self.exists:
            return f
        self.calc.start_timing('f_pp')
        lst = self.calc.el.get_property_lists(['i','s'])
        Rijn = self.calc.el.rijn
        for i,si in lst:
            for j,sj in lst:
                V = self._get_full_v(i,j)
                for n,rijn in enumerate(Rijn[:,i,j]):
                    if i==j and n==0: continue
                    dijn = np.sqrt( rijn[0]**2+rijn[1]**2+rijn[2]**2 )
                    if dijn<self.rcut:
                        f[i,:] = f[i,:] + V(dijn,der=1)*rijn/dijn
        self.calc.stop_timing('f_pp')
        return f


    def get_table(self, i, j, n=1000):
        """
        Tabulate the pair potential and return the table


        parameters:
        ===========
        i,j:    * atom indices, if integers (0,1,2,...)
                * elements, if strings ('C','H',...)
        n:      Number of grid points
        """
        if isinstance(i,int):
            pp = self.atomv[ij2s(i,j)]
        else:
            pp = self.elemv[i+j]

        Rmax = pp(None)

        x   = np.linspace(0.0, Rmax, n)
        y   = pp(x, der=0)
        dy  = pp(x, der=1)

        return x, y, dy

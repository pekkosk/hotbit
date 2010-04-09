# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

import numpy as nu
from ase.units import Bohr,Hartree
from box import mix
#from box import Atoms
from box.interpolation import Function
import os
import ase
from weakref import proxy
#import cgitb; cgitb.enable()
find=mix.find_value
vec=nu.array
norm=nu.linalg.norm
sqrt=nu.sqrt

class Repulsion:
    def __init__(self,calc):
        self.calc=proxy(calc)

        # read repulsions
        self.vrep={}
        present=calc.el.get_present()
        self.files=self.calc.ia.get_files()
        self.rmax=0.0
        for si in present:
            for sj in present:
                try:
                    table_ij=mix.find_value(self.files[si+sj],'%s_%s_table' %(si,sj),fmt='matrix')
                    table_ji=mix.find_value(self.files[sj+si],'%s_%s_table' %(sj,si),fmt='matrix')
                except KeyError:
                    raise KeyError('Interaction file for %s-%s or %s-%s not found.' %(si,sj,sj,si))
                self.vrep[si+sj]=RepulsivePotential(self.files[si+sj])
                self.rmax = max( self.rmax, self.vrep[si+sj].get_r_cut() )
        self.N=calc.el.get_N()

    def __del__(self):
        pass

    def greetings(self):
        """ Return the repulsion documentations. """
        txt='Repulsions:\n'
        shown=[]
        for ia in self.files:
            file=self.files[ia]
            if file in shown: continue
            txt+='  %s in %s\n' %(ia,file)
            doc=mix.find_value(file,'repulsion_comment',fmt='strings',default=['no repulsion doc'])
            for line in doc:
                txt+='    *'+line.lstrip()+'\n'
            shown.append(file)
        return txt
    
    
    def get_repulsion_distances(self,ela,elb,r_cut):
        """ 
        Return distances below r_cut for given elements.

        Intended use for repulsion fitting. Go through all
        periodic images, if any (of course).

        parameters:
        ===========
        ela:        chemical symbol a
        elb:        chemical symbol b
        r_cut:      repulsion cutoff (Angstrom)

        return:
        =======
        d:          numpy array of distances (in Angstroms) below r_cut
        """
        lst = self.calc.el.get_property_lists(['i','s'])
        Rijn = self.calc.el.rijn
        r_cut /= Bohr
        
        dlist=[]
        for i,si in lst:
            for j,sj in lst[i:]:
                if not ((si==ela and sj==elb) or (sj==ela and si==elb)) :
                    continue
                for n,rijn in enumerate(Rijn[:,i,j]):
                    if i==j and n==0: continue
                    d = nu.sqrt( rijn[0]**2+rijn[1]**2+rijn[2]**2 )
                    if d<=r_cut:
                        dlist.append(d)
        return nu.array(dlist)*Bohr
                    
    

    def get_repulsive_energy(self):
        """ Return the repulsive energy (in eV). """
        self.calc.start_timing('e_rep')
        
        # TODO: be more selective with n (for efficiency)
        lst=self.calc.el.get_property_lists(['i','s'])
        Rijn = self.calc.el.rijn
        erep=0.0
        for i,si in lst: 
            for j,sj in lst[i:]:
                for n,rijn in enumerate(Rijn[:,i,j]): 
                    if i==j and n==0: continue
                    d = nu.sqrt( rijn[0]**2+rijn[1]**2+rijn[2]**2 )
                    if d>self.rmax: continue
                    # TODO: remove following assert
                    assert(d>1E-10)
                    if i==j:
                        erep+=0.5*self.vrep[si+sj](d)
                    else:
                        erep+=self.vrep[si+sj](d)
        self.calc.stop_timing('e_rep') 
        return erep*Hartree
    
    
    def get_pair_repulsive_energy(self,i,j):
        """
        Return the repulsive energy of given atom pair.
        
        parameters:
        ===========
        i,j:     atom indices
        """
        si = self.calc.el.symbols[i]
        sj = self.calc.el.symbols[j]                        
        Rijn = self.calc.el.rijn
        erep=0.0
        for n,rijn in enumerate(Rijn[:,i,j]): 
            if i==j and n==0: continue
            d = nu.sqrt( rijn[0]**2+rijn[1]**2+rijn[2]**2 )
            if d>self.rmax: continue
            if i==j:
                erep += 0.5*self.vrep[si+sj](d)
            else:
                erep += self.vrep[si+sj](d)      
        return erep
    

    def get_repulsive_forces(self):
        """ 
        Return repulsive forces. 
        
        F_i = sum_(j,n) V'(rijn) rijn/dijn, with rijn = r_j^n -r_i and dijn=|rijn|
        """
        self.calc.start_timing('f_rep')
        f=nu.zeros((self.N,3))
        lst = self.calc.el.get_property_lists(['i','s'])
        Rijn = self.calc.el.rijn
        for i,si in lst:
            for j,sj in lst:
                V = self.vrep[si+sj]
                for n,rijn in enumerate(Rijn[:,i,j]):
                    if i==j and n==0: continue
                    dijn = sqrt( rijn[0]**2+rijn[1]**2+rijn[2]**2 )
                    if dijn<self.rmax:
                        f[i,:] = f[i,:] + V(dijn,der=1)*rijn/dijn
        self.calc.stop_timing('f_rep')
        return f


class RepulsivePotential:
    def __init__(self,repulsion):
        if type(repulsion)==type(''):
            self.read_repulsion(repulsion)
        else:
            self.v=repulsion
        self.range=None
            
    def get_r_cut(self):
        """ Return the cut-off of the potential. """
        return self.r_cut

    def __call__(self,r,der=0):
        """ Return V_rep(r) or V_rep'(r) """
        if r>self.r_cut:
            return 0.0
        else:
            return self.v(r,der=der)

    def read_repulsion(self,file):
        """ Read the repulsive potential from par-file. """
        try:
            v=mix.find_value(file,'repulsion',fmt='matrix')
        except:
            v=nu.array([[0,0],[1,0],[2,0],[3,0]])
        self.v=Function('spline',v[:,0],v[:,1],k=3)
        self.r_cut=v[-1,0]

    def plot(self):
        """ Plot vrep and derivative. """
        import pylab as pl
        rmin=0.5*self.r_cut
        r=nu.linspace(rmin,self.r_cut)
        v=vec([self(x,der=0) for x in r])
        vp=vec([self(x,der=1) for x in r])

        # Vrep
        pl.subplots_adjust(wspace=0.25)
        pl.subplot(1,2,1)
        pl.ylabel('$V_{rep}(r) (Ha)$')
        pl.xlabel('$r (Bohr)$')
        pl.axvline(x=self.r_cut,c='r',ls=':')
        pl.plot(r,v)
        pl.ylim(ymin=0,ymax=self(rmin))
        print 'ymax',self(rmin)
        pl.xlim(xmin=rmin)

        # Vrep'
        pl.subplot(1,2,2)
        pl.ylabel('$dV_{rep}(r)/dr (Ha/Bohr)$')
        pl.xlabel('$r (Bohr)$')
        pl.plot(r,vp,label='$V_{rep}(r)$')
        pl.axvline(x=self.r_cut,c='r',ls=':')
        pl.xlim(xmin=rmin)
        pl.ylim(ymin=self(rmin,der=1))
        pl.legend()
        pl.show()





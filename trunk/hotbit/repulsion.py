# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

import numpy as nu
from ase.units import Bohr,Hartree
from box import mix
from box import Atoms
from box.interpolation import Function
import os
import ase
from weakref import proxy
#import cgitb; cgitb.enable()
find=mix.find_value
vec=nu.array

class Repulsion:
    def __init__(self,calc):
        self.calc=proxy(calc)
        self.timer=calc.timer

        # read repulsions
        self.vrep={}
        present=calc.el.get_present()
        self.files=self.calc.ia.get_files()
        for si in present:
            for sj in present:
                try:
                    table_ij=mix.find_value(self.files[si+sj],'%s_%s_table' %(si,sj),fmt='matrix')
                    table_ji=mix.find_value(self.files[sj+si],'%s_%s_table' %(sj,si),fmt='matrix')
                except KeyError:
                    raise KeyError('Interaction file for %s-%s or %s-%s not found.' %(si,sj,sj,si))
                self.vrep[si+sj]=RepulsivePotential(self.files[si+sj])

    def __del__(self):
        print "Repulsion deleted"

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

    def get_repulsive_energy(self):
        """ Calculate the repulsive energy with included element pairs. """
        erep=0.0
        for i,j,si,sj,dist in self.calc.el.get_ia_atom_pairs(['i','j','si','sj','dist']):
            if i==j:
                continue
            erep+=self.vrep[si+sj](dist)
        return erep*Hartree

    def get_repulsive_forces(self):
        """ Calculate the forces due to repulsive potentials for element pairs. """
        self.timer.start('f_rep')
        f=nu.zeros((len(self.calc.el),3))
        for i,j,si,sj,dist,rhat in self.calc.el.get_ia_atom_pairs(['i','j','si','sj','dist','rhat']):
            if i==j:
                continue
            frep=self.vrep[si+sj](dist,der=1)*rhat
            f[i,:]=f[i,:]+frep
            f[j,:]=f[j,:]-frep
        self.timer.stop('f_rep')
        return f




class RepulsivePotential:
    def __init__(self,repulsion):
        if type(repulsion)==type(''):
            self.read_repulsion(repulsion)
        else:
            self.v=repulsion

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





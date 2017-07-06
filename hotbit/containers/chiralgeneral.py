from __future__ import division

import numpy as np
from box.mix import phival
from math import sin,cos
from weakref import proxy
import warnings

class ChiralGeneral:

    def __init__(self,atoms,type):
        '''
        Class for general chiral boundary conditions with two operations.

        @param: atoms    hotbit.Atoms instance
        @param: type     Should equal to "ChiralGeneral"

        More documentation for the methods can be found from hotbit.Atoms -class.
        '''
        self.type='ChiralGeneral'
        assert type==self.type
        self.atoms = proxy(atoms)
        self.par = {'trans1':(1,0),'angle1':(2,0),'trans2':(0,1),'angle2':(1,1),'physical':(2,1)}
        self.atoms.set_pbc((True,True,False))

    def get_type(self):
        return self.type

    def __repr__(self):
        trans1, angle1, trans2, angle2, physical = self.get('trans1'),self.get('angle1'),self.get('trans2'),self.get('angle2'),self.get('physical')
        x='ChiralGeneral: trans1=%.4f, angle1=%.4f\n' %(trans1,angle1)
        x+='               trans2=%.4f, angle2=%.4f\n' %(trans2,angle2)
        if physical:
            x+='physical), '
        else:
            x+='not physical), '
        return x

    def get_table(self):
        M = np.Inf
        return [{'M':M},{'M':M},{'M':1}]

    def get(self,key):
        """
        Get container parameters

        key: 'tran1','angle1','trans2','angle2','physical'
        """
        x = self.atoms.get_cell()[self.par[key]]
        if key in ['trans1','angle1','trans2','angle2']:
            return x
        else:
            return bool(np.round(x))

    def _set(self,**kwargs):
        assert len(kwargs)==1
        for key in kwargs:
            cell = self.atoms.get_cell()
            cell[self.par[key]] = kwargs[key]
            self.atoms.set_cell(cell)


    def set(self, trans1=None, angle1=None, trans2=None, angle2=None, physical=False, scale_atoms=False, container=None):
        """

        parameters:
        ===========
        trans1   translation in chiral operation 1
        angle1   twist angle in chiral operation 1
        trans2   translation in chiral operation 2
        angle2   twist angle in chiral operation 2
        physical if angle is small, it does not be
                         exactly 2*pi/integer, i.e. situation has no physical meaning
                         (use for calculating stuff continuously)
        scale_atoms Scale atoms according to changes in parameters
        """
        if container!=None:
            assert angle==None and height==None and M==None and twist==None
            self.set(trans1=container.get('trans1'),angle1=container.get('angle1'),\
                     trans2=container.get('trans2'),angle2=container.get('angle2'),\
                     physical=container.get('physical'))

        if trans1!=None: self._set(trans1=trans1)
        if angle1!=None: self._set(angle1=angle1)
        if trans2!=None: self._set(trans2=trans2)
        if angle2!=None: self._set(angle2=angle2)
        if physical:
            raise NotImplementedError('Physical=True not implemented for ChiralGeneral.')
        if physical!=None: self._set(physical=physical)
        if scale_atoms:
            raise NotImplementedError('Scaling of atom positions not implemented in ChiralGeneral.')

    def __eq__(self,other):
        return self.atoms == other.atoms


    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        M = np.Inf
        ranges = np.array([[-M,M],[-M,M],[0,0]])
        return ranges


    def transform(self,r,n):
        """ Rotate around z r by (n2*angle+n0*twist) and translate by n0*height. """
        R = self.rotation(n)
        trans = np.zeros((3))
        trans = n[0]*np.array([0,0,self.get('trans1')]) + n[1]*np.array([0,0,self.get('trans2')])
        return np.dot(R,r) + np.array(trans)

    def rotation(self,n,angles=False):
        """ Active rotation matrix of given angle wrt. z-axis."""
        angle = n[0]*self.get('angle1') + n[1]*self.get('angle2')
        R = np.array([[cos(angle),-sin(angle),0],[sin(angle),cos(angle),0],[0,0,1]])
        if angles:
            return 0.,0.,angle
        else:
            return R

from __future__ import division

import numpy as np
from box.mix import phival
from math import sin,cos
from weakref import proxy
import warnings

class Wedge:

    def __init__(self,atoms,type):
        '''
        Class for wedge boundary conditions.

           ^ y-axis
           |    /
           |   /
           |  /
           | /  angle
           +----------> x-axis

        @param: atoms    hotbit.Atoms instance
        @param: type     Should equal to "Wedge"

        More documentation for the methods can be found from hotbit.Atoms -class.
        '''
        self.type='Wedge'
        assert type==self.type
        self.atoms = proxy(atoms)
        self.par = {'height':(1,0),'angle':(2,0),'pbcz':(0,1),'physical':(1,1)}
        self.atoms.set_pbc((True,False,atoms.get_pbc()[2]))
        #self._set_table()

    def get_type(self):
        return self.type

    def __repr__(self):
        angle, height, pbcz, physical = self.get('angle'), self.get('height'), self.get('pbcz'), self.get('physical')
        x='Wedge: angle=%.4f (2*pi/%.2f, ' %(angle,2*np.pi/angle)
        if physical:
            x+='physical), '
        else:
            x+='not physical), '
        x+='height=%.4f Ang ' %height
        if pbcz:
            x+='(z:pbc)'
        else:
            x+='(z:no pbc)'
        return x

    def get_table(self):
        M = int( round(2*np.pi/self.get('angle')) )
        if self.get('pbcz'):
            return [{'M':M},{'M':1},{'M':np.Inf}]
        else:
            return [{'M':M},{'M':1},{'M':1}]


    def get(self,key):
        """
        Get container parameters

        key: 'angle','height','pbcz','physical'
        """
        if key=='pbcz':
            return self.atoms.get_pbc()[2]
        else:
            x = self.atoms.get_cell()[self.par[key]]
            if key in ['angle','height']:
                return x
            else:
                return bool(np.round(x))

    def _set(self,**kwargs):
        assert len(kwargs)==1
        if 'pbcz' in kwargs:
            self.atoms.set_pbc( (True,False,kwargs['pbcz']) )
        else:
            for key in kwargs:
                cell = self.atoms.get_cell()
                cell[self.par[key]] = kwargs[key]
                self.atoms.set_cell(cell)


    def set(self, angle=None, height=None, M=None, physical=True, pbcz=None, scale_atoms=False, container=None):
        """ Only height can be reset, not angle.

        parameters:
        ===========
        angle    angle (in radians) of the wedge (and M=None)
        height   Height of the primitive cell in z-direction
        M        set angle to 2*pi/M (and angle=None)
        physical (only if M=None) if angle is small, it does not be
                         exactly 2*pi/integer, i.e. situation has no physical meaning
                         (use for calculating stuff continuously)
        pbcz     True if wedge is periodic in z-direction
        scale_atoms Scale atoms according to changes in parameters. When changing angle, scale
                    also radial distances. (Radii are inversely proportional to angle.)
        """
        if container!=None:
            assert angle==None and height==None and M==None and pbcz==None
            self.set(angle=container.get('angle'),height=container.get('height'),\
                     physical=container.get('physical'), pbcz=container.atoms.get_pbc()[2])

        if angle!=None or M!=None:
            #assert not scale_atoms
            assert not (angle!=None and M!=None)
            old_angle = self.get('angle')
            if M != None:
                assert isinstance(M,int)
                self._set(angle=2*np.pi/M)
            elif angle != None:
                M = int( round(2*np.pi/angle) )
                self._set(angle=angle)

            # check parameters
            self._set( physical=float(physical) )
            if self.get('angle')<1E-6:
                raise Warning('Too small angle (%f) may bring numerical problems.' %self.get('angle'))
            if self.get('angle')>np.pi:
                raise AssertionError('angle>pi')
            if np.abs(M-2*np.pi/self.get('angle'))>1E-12 and self.get('physical'):
                raise AssertionError('angle not physical: angle != 2*pi/M')
            if not self.get('physical') and M<20:
                warnings.warn('Quite large, non-physical angle 2*pi/%.4f.' %(2*np.pi/self.get('angle')) )

            if scale_atoms:
                if abs(old_angle)<1E-10:
                    raise ValueError('Atoms cannot be scaled; old wedge angle too small.')
                newr = []
                for r in self.atoms.get_positions():
                    x,y = r[0],r[1]
                    rad = np.sqrt( x**2+y**2 )
                    newphi = phival(x,y)*(self.get('angle')/old_angle)
                    rad2 = rad * old_angle/self.get('angle')
                    newr.append( [rad2*np.cos(newphi),rad2*np.sin(newphi),r[2]] )
                self.atoms.set_positions(newr)

        if height!=None:
            if scale_atoms:
                r = self.atoms.get_positions()
                r[:,2] = r[:,2] * height/self.get('height')
                self.atoms.set_positions(r)
            self._set(height=height)

        if pbcz!=None:
            self._set(pbcz=float(pbcz))

        #self._set_table()


    def __eq__(self,other):
        return self.atoms == other.atoms


    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        M = int( round(2*np.pi/self.get('angle')) )
        i = M//2
        zi = 0
        if self.get('pbcz'):
            zi = np.Inf
        if np.mod(M,2)==1:
            ranges = np.array([[-i,i],[0,0],[-zi,zi]])
        else:
            ranges = np.array([[-i+1,i],[0,0],[-zi,zi]])
        return ranges


    def transform(self,r,n):
        """ Rotate r by n2*angle. """
        R = self.rotation(n)
        trans = np.zeros((3))
        if self.get('pbcz'):
            trans = n[2]*np.array([0,0,self.get('height')])
        return np.dot(R,r) + np.array(trans)

    def rotation(self,n,angles=False):
        """ Active rotation matrix of given angle wrt. z-axis."""
        angle = n[0]*self.get('angle')
        if angles:
            return 0.,0.,angle
        else:
            R = np.array([[cos(angle),-sin(angle),0],[sin(angle),cos(angle),0],[0,0,1]])
            return R

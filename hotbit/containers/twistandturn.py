import numpy as np
from math import sin,cos
from weakref import proxy
from box.mix import rotation_matrix, rotation_from_matrix, phival

class TwistAndTurn:
    def __init__(self,atoms,type):
        '''
        Class for bent chiral boundary conditions.

        atoms:    hotbit.Atoms -instance
        type:     Should equal to "TwistAndTurn"

        Periodic direction is the third symmetry operation.

        More documentation for the methods can be
        found from hotbit.Atoms -class.
        '''
        self.type='TwistAndTurn'
        assert type == self.type
        self.atoms = proxy(atoms)
        self.R = None
        self.baxis = np.array([0,0,1]) # bending axis
        self.par = {'bend_angle':(0,1),'twist_angle':(0,2),'R':(1,2)}
        self._set_table()
        self.atoms.set_pbc( (False,False,True) )

    def _set_table(self):
        self.table = [{'M':1},{'M':1},{'M':np.inf}]

    def get_type(self):
        return self.type

    def __repr__(self):
        a1,a2,R = [self.get(k) for k in ['bend_angle','twist_angle','R']]
        x='TwistAndTurn: bend_angle=%.4f, twist_angle=%.4f, R=%.4f' %(a1,a2,R)
        return x

    def get_table(self):
        return [{'M':1},{'M':1},{'M':np.Inf}]

    def get(self, key):
        """
        Return container parameters.

        parameters:
        ===========
        key:    'bend_angle','twist_angle','R'
        """
        cell = self.atoms.get_cell()
        if key in ['bend_angle','twist_angle']:
            x = cell[self.par[key]]
        elif key=='R':
            x = cell[self.par[key]]
            # adapt R to atoms' mean radius
            if np.abs(x)<1E-16:
                r = self.atoms.get_positions()
                x = np.sqrt(r[:,0]**2+r[:,1]**2).mean()
        return x

    def _set(self,**kwargs):
        """
        Save parameters to cell matrix
        """
        assert len(kwargs)==1
        cell = self.atoms.get_cell()
        for key in kwargs:
            cell[self.par[key]] = kwargs[key]
        self.atoms.set_cell(cell)


    def set(self,bend_angle=None,twist_angle=None,R=None,mode=None,container=None,scale_atoms=False):
        """
        Reset angles, R or axes

        @param: bend_angle   bending angle
        @param: twist_angle  twisting angle
        @param: R            radius of curvature for neutral line
                             (if R==0, use R based on mean distance of atoms)
        @param: scale_atoms  scale atoms according to container modifications
        """
        if container is not None:
            # copy container
            assert bend_angle is None
            assert twist_angle is None
            assert R is None
            assert mode is None
            self.set( angle1=container.bend_angle,
                      angle2=container.twist_angle,
                      R=container.R,
                      mode=container.mode)
        if bend_angle is not None:
            if scale_atoms:
                ratio = bend_angle/self.get('bend_angle')
                r = self.atoms.get_positions()
                R0 = np.sqrt( r[:,0]**2+r[:,1]**2 )
                angle = np.array([np.arctan2(y,x) for x,y in zip(r[:,0],r[:,1])] )
                angle2 = ratio*angle
                r2 = np.array( [R0*np.cos(angle2),R0*np.sin(angle2),r[:,2]] ).transpose()
                self.atoms.set_positions(r2)
            self._set(bend_angle = bend_angle)
        if twist_angle is not None:
            if scale_atoms:
                da = twist_angle - self.get('twist_angle')
                r = self.atoms.get_positions()
                angle = np.array([np.arctan2(y,x) for x,y in zip(r[:,0],r[:,1])] )
                x = angle/self.get('bend_angle')
                r2 = np.array( [self.transform_arbitrary(r[i],0,da*x[i]) for i in range(len(self.atoms))] )
                self.atoms.set_positions(r2)
            self._set(twist_angle = twist_angle)
        if R is not None:
            if scale_atoms:
                dR = R-self.get('R')
                r = self.atoms.get_positions()
                R0 = np.sqrt( r[:,0]**2+r[:,1]**2 )
                angle = np.array([np.arctan2(y,x) for x,y in zip(r[:,0],r[:,1])] )
                r[:,0] = (R0+dR)*np.cos(angle)
                r[:,1] = (R0+dR)*np.sin(angle)
                self.atoms.set_positions(r)
            self._set(R = R)
        if mode is not None:
            self._set(mode = mode)

        self.atoms.set_pbc((False,False,True))


    def __eq__(self,other):
        return self.atoms == other.atoms


    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        ranges = np.array( [[0,0],[0,0],[-np.Inf,np.Inf]] )
        return ranges

    def transform_arbitrary(self,r,bend,twist):
        """ The basic symmetry transformation with arbitrary angles b (bend angle) and t (twist angle) """
        r = np.array(r)
        R = self.get('R')
        Rp = np.array([r[0],r[1],0.])     # vector on xy-plane
        Rv = R*Rp/np.linalg.norm(Rp)      # vector to 'center' of tube
        rho = r-Rv                        # vector from tube center to r
        t = np.cross(np.array([0,0,1]),r) # tangential vector

        Rtwist = rotation_matrix(t,twist)
        Rbend  = rotation_matrix(self.baxis,bend)
        return np.dot( Rbend,Rv+np.dot(Rtwist,rho) )

    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        if n[2]==0:
            return r.copy()

        bend, twist = [self.get(k) for k in ['bend_angle','twist_angle']]
        return self.transform_arbitrary(r,n[2]*bend,n[2]*twist)

    def rotation(self,n,angles=False):
        """
        Rotate around two axes:
        first around 'unit cell axis' for chirality
        and then around 'torus axis' for twisting..
        Returns the rotation matrix.

        If angles==True, return (theta,phi,angle), where (theta,phi)
        gives the rotation axis and 'angle' the rotation angle.

        Approximate tangential vector by a mean tangential vector
        (-sin(angle/2),cos(angle/2),0)
        """
        if n[2]==0:
            return np.eye(3)
        bend,twist = [self.get(k) for k in ['bend_angle','twist_angle']]
        rot1 = rotation_matrix(self.baxis,n[2]*bend)
        a = self.get('bend_angle')
        taxis = np.array([-sin(a/2),cos(a/2),0])
        rot2 = rotation_matrix(taxis,n[2]*twist)
        R = np.dot(rot2,rot1)
        if angles:
            angle, dir = rotation_from_matrix(R)
            theta = np.arccos(dir[2])
            if np.abs(theta)<1E-12 or np.abs(theta-np.pi)<1E-12:
                phi = 0.
            else:
                cosp, sinp = dir[0]/np.sin(theta), dir[1]/np.sin(theta)
                phi = phival(cosp,sinp)
            return (theta,phi,angle)
        else:
            return R

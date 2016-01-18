import numpy as np
from box.mix import phival
from math import sin,cos
from weakref import proxy
from box.mix import rotation_matrix

class TwistAndTurn:
    def __init__(self,atoms,type):
        '''
        Class for bent chiral boundary conditions.
        
        @param: atoms    hotbit.Atoms -instance
        @param: type     Should equal to "TwistAndTurn"

        More documentation for the methods can be
        found from hotbit.Atoms -class.
        '''
        self.type='TwistAndTurn'
        assert type == self.type
        self.atoms = proxy(atoms)
        self.angle1 = None
        self.angle2 = None
        self.R = None
        self.n1=([0,0,1])
        self.n2=None
        self.mode=None
        self.par = {'angle1':(0,1),
                    'angle2':(0,2),
                    'R':(1,2),
                    'n2':(1,slice(0,2)),
                    'mode':(2,0)}
        self._set_table()
        self.atoms.set_pbc( (False,False,True) )

    def _set_table(self):
        self.table = [{'M':1},{'M':1},{'M':np.inf}]

    def get_type(self):
        return self.type

    def __repr__(self):
        a1,a2,R = [self.get(k) for k in ['angle1',
                                         'angle2',
                                         'R']]
        x='TwistAndTurn: angle1=%.4f, angle2=%.4f, R=%.4f' %(a1,
                                                           a2,
                                                           R)
        return x            

    def get_table(self):
        return [{'M':1},{'M':1},{'M':np.Inf}]

    def get(self, key):
        """
        Return container parameters.
    
        parameters:
        ===========
        key:    'angle1','angle2','R','n2','mode'
        """
        cell = self.atoms.get_cell()
        x = cell[self.par[key]]
        if key == 'n2':
            return np.concatenate((x,[0]))
        else:
            return x    

    def _set(self,**kwargs):
        """
        Save parameters to cell matrix
        """
        assert len(kwargs)==1
        for key in kwargs:
            cell = self.atoms.get_cell()
            value = kwargs[key]
            if key == 'n2':
                value = value[0:2]
            cell[self.par[key]] = value
            self.atoms.set_cell(cell)
   

    def set(self,
            angle1=None,
            angle2=None,
            R=None,
            n2=None,
            mode=None,
            container=None):
        """
        Reset angles, R or axes
        
        @param: angle1   Angle for chirality
        @param: angle2   Angle for bending
        @param: R        Radius of curvature
        @param: n2       Axis for bending
        @param: mode     Mode for transformation
        """
        if container!=None:
            # copy container
            assert angle1==None
            assert angle2==None
            assert R==None
            assert n2==None
            assert mode == None
            self.set( angle1=container.angle1,
                      angle2=container.angle2,
                      R=container.R,
                      n2=container.n2,
                      mode=container.mode)
        if angle1!=None:
            self._set(angle1 = angle1)
        if angle2!=None:
            self._set(angle2 = angle2)
        if R!=None:
            self._set(R = R)
        if n2!=None:
            self._set(n2 = n2)
        if mode!=None:
            self._set(mode = mode)

        if self.get('n2')!=None:
            self._set(n2 = self.get('n2')/
                      np.linalg.norm(self.get('n2')))
        self.atoms.set_pbc((False,False,True))
        

    def __eq__(self,other):
        return self.atoms == other.atoms


    def get_symmetry_operation_ranges(self):
        """ Return ranges for symmetry operations. """
        ranges = np.array( [[0,0],[0,0],[-np.Inf,np.Inf]] )
        return ranges

    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        if n[2]==0:
            return r.copy()

        a1, a2, R, n2, mode = [self.get(k) for k in ['angle1',
                                                     'angle2',
                                                     'R',
                                                     'n2',
                                                     'mode']]
        
        n1 = self.n1
        
        if mode == 1:
            orig_chir = np.cross(n2,n1)
            orig_chir = orig_chir/np.linalg.norm(orig_chir)
            A = np.dot(R,orig_chir)

            rot1 = rotation_matrix(n1,n[2]*a1)
            rot2 = rotation_matrix(n2,-n[2]*a2)

            Rot = np.dot( rot2,rot1 )

            rp = r.copy()
            rp = np.dot( Rot,rp-A ) + np.dot( rot2,A )
            return rp

        elif mode == 2:
            x = r[0]
            y = r[1]
            z = r[2]

            alpha1 = phival(x,z)
            Rr = np.sqrt(x**2 + z**2)
            R1 = R*np.array([np.cos(alpha1),
                             0.0,
                             np.sin(alpha1)])
            rho = np.linalg.norm(r-R1)
            beta1 = phival(Rr - R, y)

            alpha2 = a2*n[2] + alpha1
            beta2 = a1*n[2] + beta1

            Rrp = R + rho*np.cos(beta2)

            xp = Rrp*np.cos(alpha2)
            yp = rho*np.sin(beta2)
            zp = Rrp*np.sin(alpha2)
            
            rp = np.array([xp,yp,zp])
            return rp

    def rotation(self,n):
        """
        Rotate around two axes:
        first around 'unit cell axis' for chirality
        and then around 'torus axis' for bending.
        Returns the rotation matrix.
        """
        if n[2]==0:
            return np.eye(3)
        a1,a2,n2 = [self.get(k) for k in ['angle1',
                                          'angle2',
                                          'n2']]
        n1 = self.n1

        rot1 = rotation_matrix(n1,n[2]*a1)
        rot2 = rotation_matrix(n2,-n[2]*a2)
       
        return np.dot(rot2,rot1)

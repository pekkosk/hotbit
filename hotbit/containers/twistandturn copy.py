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
        self.R = None
        self.baxis = np.array([0,0,1]) # bending axis
        self.taxis = np.array([0,1,0]) # (primary) twisting axis
        self.par = {'bend_angle':(0,1),'twist_angle':(0,2),'R':(1,2),'mode':(2,0)}
        self._set_table()
        self.atoms.set_pbc( (False,False,True) )

    def _set_table(self):
        self.table = [{'M':1},{'M':1},{'M':np.inf}]

    def get_type(self):
        return self.type

    def __repr__(self):
        a1,a2,R,m = [self.get(k) for k in ['bend_angle','twist_angle','R','mode']]
        x='TwistAndTurn: bend_angle=%.4f, twist_angle=%.4f, R=%.4f, mode=%i' %(a1,a2,R,int(m))
        return x            

    def get_table(self):
        return [{'M':1},{'M':1},{'M':np.Inf}]

    def get(self, key):
        """
        Return container parameters.
    
        parameters:
        ===========
        key:    'bend_angle','twist_angle','R','mode'
        """
        cell = self.atoms.get_cell()
        if key in ['bend_angle','twist_angle','R']:
            x = cell[self.par[key]]
        elif key=='mode':
            x = int( round(cell[self.par[key]]) )
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
        @param: mode     	 mode for transformation (different modes of approximation)
        @param: scale_atoms  scale atoms according to container modifications
        """
        if container!=None:
            # copy container
            assert bend_angle==None
            assert twist_angle==None
            assert R==None
            assert mode == None
            self.set( angle1=container.bend_angle,
                      angle2=container.twist_angle,
                      R=container.R,
                      mode=container.mode)
        if bend_angle!=None:
            if scale_atoms:
                ratio = bend_angle/self.get('bend_angle')
                r = self.atoms.get_positions()
                R = np.sqrt( sum(r[:,0]**2+r[:,1]**2) )
                alpha1 = np.array([phival(x1,y1) for x1,y1 in zip(r[:,0],r[:,1])] )
                alpha = ratio*alpha1
                r2 = np.array( [R*np.cos(alpha),R*np.sin(alpha),r[:,2]] ).transpose()
                self.atoms.set_positions(r2)
            self._set(bend_angle = bend_angle)
        if twist_angle!=None:
            if scale_atoms:
                da = twist_angle - self.get('twist_angle')
                R = self.get('R')
                r = self.atoms.get_positions()
                alpha1 = np.array([phival(x1,y1) for x1,y1 in zip(r[:,0],r[:,1])] )
                Rr = np.sqrt(r[:,0]**2+r[:,1]**2)
                R1 = R * np.array( [np.cos(alpha1),np.sin(alpha1),np.zeros_like(alpha1)] ).transpose()
                rho = np.sqrt( np.sum((r-R1)**2,axis=1) ) 
                beta1 = np.array( [phival(a,b) for a,b in zip(Rr-R,-r[:,2])] )

                # scale twisting angle depending on angle between x and z (alpha1)
                x = alpha1/self.get('bend_angle')
                beta2  = beta1  + x*da 
                Rrp = R + rho*np.cos(beta2)
                r2 = np.array( [Rrp*np.cos(alpha1),Rrp*np.sin(alpha1),-rho*np.sin(beta2)] ).transpose()
                self.atoms.set_positions(r2)
            self._set(twist_angle = twist_angle)
        if R!=None:
            if scale_atoms:
                dR = R-self.get('R')
                r = self.atoms.get_positions()
                rad = np.sqrt( r[:,0]**2+r[:,1]**2 )
                phi = np.array([phival(x,y) for x,y in zip(r[:,0],r[:,1])] )
                r[:,0] = (rad+dR)*np.cos(phi)
                r[:,1] = (rad+dR)*np.sin(phi)
                self.atoms.set_positions(r)
            self._set(R = R)
        if mode!=None:
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
        Rp = np.array([r[0],r[1],0.])
        Rv = R*Rp/np.linalg.norm(Rp)
        rho = r-Rv
        t = np.cross(np.array([0,0,1]),r)
        Rtwist = rotation_matrix(t,twist)
        Rbend = rotation_matrix(self.baxis,bend)
        
        return np.dot( Rbend(Rv+np.dot(Rtwist,rho)),r )

    def transform(self,r,n):
        """ Symmetry transformation n for position r. """
        if n[2]==0:
            return r.copy()

        a1, a2, R, mode = [self.get(k) for k in ['bend_angle','twist_angle','R','mode']]
        
        n1 = self.n1
        
        if mode == 2:
            orig_chir = np.cross(n2,n1)
            orig_chir = orig_chir/np.linalg.norm(orig_chir)
            A = np.dot(R,orig_chir)

            rot1 = rotation_matrix(n1,n[2]*a1)
            rot2 = rotation_matrix(n2,-n[2]*a2)

            Rot = np.dot( rot2,rot1 )

            rp = r.copy()
            rp = np.dot( Rot,rp-A ) + np.dot( rot2,A )
            return rp

        elif mode == 1:
            x, y, z = r
            x,y,z = x,z,-y
            alpha1 = phival(x,y)
            Rr = np.sqrt(x**2+y**2)
            R1 = R * np.array([np.cos(alpha1),np.sin(alpha1),np.zeros_like(alpha1)])
            rho = np.linalg.norm(r-R1)
            beta1 = phival(Rr-R,z)

            alpha2 = alpha1 + n[2]*a1 
            beta2  = beta1  + n[2]*a2 

            Rrp = R + rho*np.cos(beta2)
            return np.array( [Rrp*np.cos(alpha2),Rrp*np.sin(alpha2),-rho*np.sin(beta2)] )
        else:
            raise AssertionError('TwistAndTurn mode is not set')

    def rotation(self,n):
        """
        Rotate around two axes:
        first around 'unit cell axis' for chirality
        and then around 'torus axis' for twisting..
        Returns the rotation matrix.
        """
        if n[2]==0:
            return np.eye(3)
        a1,a2 = [self.get(k) for k in ['bend_angle','twist_angle']]
        rot1 = rotation_matrix(self.baxis,n[2]*a1)
        rot2 = rotation_matrix(self.taxis,n[2]*a2)
        return np.dot(rot2,rot1)

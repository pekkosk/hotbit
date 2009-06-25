from ase import Atoms as ase_Atoms
import numpy as nu
from box.mix import  phival
from math import sin,cos
from box import mix


class WedgeAtoms(ase_Atoms):     
    
    def __init__(self, *args, **kwargs):
        ase_Atoms.__init__(self, *args, **kwargs)
        self.pbc=nu.array([True,False,False],bool)
               
               
    def set(self,omega,physical=True):
        # TODO: disable set_pbc and set_cell
        self.omega=omega
        M1 = round(2*nu.pi/omega)
        if omega<1E-6:
            raise Warning('Too small omega may bring rounding problems.')
        if omega>nu.pi:
            raise AssertionError('Omega>2*pi')
        if nu.abs(M1-2*nu.pi/omega)>1E-12 and physical: 
            raise AssertionError('omega not physical 2*pi/M !')
        if not physical and M1<20:
            raise AssertionError('Too large, non-physical omega.')
        self.M1=int(M1)
        i = M1/2
        if nu.mod(M1,2)==1:
            self.ranges = nu.array([[-i,i],[0,0],[0,0]],int)
        else:
            self.ranges = nu.array([[-i+1,i],[0,0],[0,0]],int)
    
        
    def multiply(self,operations):
        atoms2=None
        for n in operations:
            self._check(n)
            atomsn=self.copy()
            atomsn.set_positions( [self.transform(r,n) for r in self.get_positions()] )
            try:
                atoms2 += atomsn
            except:
                atoms2 = atomsn
        return atoms2    
    
    
    def get_ranges(self):
        '''
        Return the ranges for symmetry operations.
        '''
        return self.ranges.copy()        
    
    
    def _check(self,n):
        '''
        Check that the symmetry operation is allowed.
        
        @param n: 3-tuple for transformation
        '''
        for i in range(3):
            if not self.ranges[i,0]<=n[i]<=self.ranges[i,1]:
                raise AssertionError('Illegal symmetry operation: %i %i %i' %n)
        
    
    def transform(self,r,n):
        '''
        Transform position r with S(n)
        
        @param r: position vector
        @param n: 3-tuple for symmetry operations
        '''
        self._check(n)
        rn=nu.zeros_like(r)
        x, y = r[0], r[1]
        rad = nu.sqrt(x**2+y**2)
        phi = phival(x,y)
        rn[0] = rad * cos( phi+n[0]*self.omega )
        rn[1] = rad * sin( phi+n[0]*self.omega )
        rn[2] = r[2]
        return rn
        
    
    def dtensor(self,r,n):
        '''
        Return the dyadic tensor 
        
                    d (R_j^n)_a
        T(jn)_ab = ------------- hat_a hat_b 
                    d (R_j)b
                    
        @param r: 
        @param n:
        '''
        self._check(n)
        rn = self.transform(r,n)
        rad = nu.sqrt(r[0]**2+r[1]**2)
        assert rad>1E-10
        T = nu.array([[r[0]*rn[0]+r[1]*rn[1],-(r[0]*rn[1]-r[1]*rn[0]),0],\
                      [r[0]*rn[1]-r[1]*rn[0],  r[0]*rn[0]+r[1]*rn[1] ,0],\
                      [       0,                      0,         rad**2]])/rad**2
        return T
        
    
    def axis_rotation(self,n):
        '''
        Return the rotation matrix for given symmetry operation.
        
        @param n: 3-tuple
        '''
        self._check(n)
        angle = n[0]*self.omega
        R = nu.array([[cos(angle),sin(angle),0],[-sin(angle),cos(angle),0],[0,0,1]])
        return R     
                       
        
    def __eq__(self,other):
        if ase_Atoms.__eq__(self,other): # and self.omega==other.omega:
            return True
        else:
            return False
    
    
    def __ne__(self,other):
        eq = self.__eq__(other)
        if eq is NotImplemented:
            return eq
        else:
            return not eq
#        
#        if isinstance(other,WedgeAtoms):
#            # More strict comparison for BravaisAtoms
#            # TODO: check additional stuff for other generalized classes
#            # (angles, torsions ...)
#            if (self.positions-other.positions<1E-13).all() and \
#               (self.pbc==other.pbc).all() and \
#               self.get_chemical_symbols()==other.get_chemical_symbols() and \
#               (self.cell==other.cell).all():
#                return True
#            else:
#                return False
#        else:
#            # other is probably normal ase.Atoms; more loose check
#            if (self.positions-other.positions<1E-13).all() and \
#               (self.pbc==other.pbc).all() and \
#               self.get_chemical_symbols()==other.get_chemical_symbols() and \
#               (self.cell==other.cell).all():
#                return True
#            else:
#                return False




class ChiralAtoms(ase_Atoms):     
    def __init__(self, *args, **kwargs):
        ase_Atoms.__init__(self, *args, **kwargs)
        self.angle = 0.0
        self.length = self.get_cell()[2,2]
        self.pbc=nu.array([False,False,True],bool)
        self.ranges=nu.array([[0,0],[0,0],[-nu.Inf,nu.Inf]])
               
               
    def set(self,angle=None,length=None):
        # TODO: disable set_pbc and set_cell
        if angle!=None:
            self.set_angle(angle)
        if length!=None:
            self.set_length(length)
        
    
    def get_ranges(self):
        '''
        Return the ranges for symmetry operations.
        '''
        return self.ranges.copy()        
    
    
    def _check(self,n):
        '''
        Check that the symmetry operation is allowed.
        
        @param n: 3-tuple for transformation
        '''
        for i in range(3):
            if not self.ranges[i,0]<=n[i]<=self.ranges[i,1]:
                raise AssertionError('Illegal symmetry operation: %i %i %i' %n)
        
    
    def transform(self,r,n):
        '''
        Transform position r with S(n)
        
        @param r: position vector
        @param n: 3-tuple for symmetry operations
        '''
        self._check(n)
        rn=nu.zeros_like(r)
        x, y = r[0], r[1]
        rad = nu.sqrt(x**2+y**2)
        phi = phival(x,y)
        rn[0] = rad * cos( phi+n[2]*self.angle )
        rn[1] = rad * sin( phi+n[2]*self.angle )
        rn[2] = r[2] + n[2]*self.length
        return rn
        
    
    def dtensor(self,r,n):
        '''
        Return the dyadic tensor 
        
                    d (R_j^n)_a
        T(jn)_ab = ------------- hat_a hat_b 
                    d (R_j)b
                    
        @param r: position vector
        @param n: symmetry operation 3-tuple
        '''
        self._check(n)
        rn = self.transform(r,n)
        rad = nu.sqrt(r[0]**2+r[1]**2)
        assert rad>1E-10
        T = nu.array([[r[0]*rn[0]+r[1]*rn[1],-(r[0]*rn[1]-r[1]*rn[0]),0],\
                      [r[0]*rn[1]-r[1]*rn[0],  r[0]*rn[0]+r[1]*rn[1] ,0],\
                      [       0,                      0,         rad**2]])/rad**2
        return T
        
    
    def axis_rotation(self,n):
        '''
        Return the rotation matrix for given symmetry operation.
        
        @param n: 3-tuple
        '''
        self._check(n)
        angle = n[2]*self.angle
        R = nu.array([[cos(angle),sin(angle),0],[-sin(angle),cos(angle),0],[0,0,1]])
        return R     
    
    
    def set_angle(self,angle,scale_atoms=True,absolute=True):
        '''
        Set the chiral angle.
        
        @param angle: the chiral angle in radians
        @param scale_atoms: scale atoms according to twisting
        @param absolute: if True, given angle is the absolute angle, otherwise
                         it is addition to present chiral angle
        '''
        if absolute:
            da = angle - self.angle
            self.angle = angle
        else:
            da = angle
            self.angle += da    
        
        if scale_atoms:
            newr = []
            for r in self.get_positions():
                x,y = r[0],r[1]
                rad = nu.sqrt( x**2+y**2 )
                newphi = mix.phival(x,y) + r[2]/self.length * da
                newr.append( [rad*nu.cos(newphi),rad*nu.sin(newphi),r[2]] )
            self.set_positions(newr)
        
        
    def set_length(self,length,scale_atoms=True):
        '''
        Set the length of the chiral system.
        
        @param length: set z-length in Angstroms
        @param scale_atoms: scale atom positions accordingly
        '''
        self.length=length
        cell = self.get_cell()
        self.set_cell([cell[0,0],cell[1,1],length],scale_atoms=scale_atoms)
        
     
    def multiply(self,operations):
        atoms2=None
        for n in operations:
            self._check(n)
            atomsn=self.copy()
            atomsn.set_positions( [self.transform(r,n) for r in self.get_positions()] )
            try:
                atoms2 += atomsn
            except:
                atoms2 = atomsn
        return atoms2
   
#        
#    def __eq__(self,other):
#        if isinstance(other,ChiralAtoms):
#            # More strict comparison for BravaisAtoms
#            # TODO: check additional stuff for other generalized classes
#            # (angles, torsions ...)
#            if (self.positions-other.positions<1E-13).all() and \
#               (self.pbc==other.pbc).all() and \
#               self.get_chemical_symbols()==other.get_chemical_symbols() and \
#               (self.cell==other.cell).all():
#                return True
#            else:
#                return False
#        else:
#            # other is probably normal ase.Atoms; more loose check
#            if (self.positions-other.positions<1E-13).all() and \
#               (self.pbc==other.pbc).all() and \
#               self.get_chemical_symbols()==other.get_chemical_symbols() and \
#               (self.cell==other.cell).all():
#                return True
#            else:
#                return False




class BravaisAtoms(ase_Atoms):     
    def __init__(self, *args, **kwargs):
        ase_Atoms.__init__(self, *args, **kwargs)
        ranges = []
        for i in range(3):
            if self.pbc[i]:
                ranges.append([-nu.Inf,nu.Inf])
            else:
                ranges.append([0,0])
        self.ranges = nu.array(ranges)
               
    
    def transform(self,r,n):
        '''
        Transform position r with S(n)
        
        @param r: position vector
        @param n: 3-tuple for symmetry operations
        '''
        # TODO: check than n belongs to allowed transformations
        rn=r.copy()
        for a in range(3):
            rn = rn + n[a]*self.cell[a,:]
        return rn
    
    
    def get_ranges(self):
        '''
        Return the ranges for symmetry operations.
        '''
        return self.ranges.copy() 
    
    
    def dtensor(self,r,n):
        '''
        Return the dyadic tensor 
        
                    d (R_j^n)_a
        T(jn)_ab = ------------- hat_a hat_b 
                    d (R_j)b
                    
        @param r: 
        @param n:
        '''
        return nu.eye(3)
    
    
    def axis_rotation(self,n):
        '''
        Return the rotation matrix for given symmetry operation.
        
        @param n: 3-tuple
        '''
        return nu.eye(3)       
        
        
#    def __eq__(self,other):
#        #print 'comp'
#        if (nu.abs(self.positions-other.positions)<1E-13).all() and \
#           (self.pbc==other.pbc).all() and \
#           self.get_chemical_symbols()==other.get_chemical_symbols() and \
#           (self.cell==other.cell).all():
#            #print 'ccccc'
#            return True
#        else:
##            print self.positions
##            print other.positions
##            print self.pbc
##            print other.pbc
##            print self.cell
##            print other.cell
#            return False
        
#        if isinstance(other,BravaisAtoms):
#            # More strict comparison for BravaisAtoms
#            # TODO: check additional stuff for other generalized classes
#            # (angles, torsions ...)
#            if (self.positions-other.positions<1E-13).all() and \
#               (self.pbc==other.pbc).all() and \
#               self.get_chemical_symbols()==other.get_chemical_symbols() and \
#               (self.cell==other.cell).all():
#                return True
#            else:
#                return False
#        else:
#            # other is probably normal ase.Atoms; more loose check
#            if (self.positions-other.positions<1E-13).all() and \
#               (self.pbc==other.pbc).all() and \
#               self.get_chemical_symbols()==other.get_chemical_symbols() and \
#               (self.cell==other.cell).all():
#                return True
#            else:
#                return False
           
        
#    def get_number_of_cells(self):

if __name__=='__main__':
    from ase import *
    atoms=Atoms('C',[(0,0,0)],cell=[(1,1,0),(0,1,0),(0,0,1)])
#    gy=NormalGeometry('C',[(0,0,0)],cell=[(1,1,0),(0,1,0),(0,0,1)])
    gy=NormalGeometry(atoms)
    print gy.vector(array([0,0,0]),n=(1,1,1),r0=[0,0,0],lst=['vec','hat','norm'])
    print gy.vector(array([0,0,0]),n=(1,1,1),lst='dtensor')
    print gy.vector(array([0,0,0]),n=(1,1,1),r0=[2,2,2],lst='vec')
        
        
        
    
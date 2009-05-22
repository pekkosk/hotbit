from ase import Atoms
from wedgeatoms import WedgeAtoms

from math import pi, cos, sin

class BendAtoms(object):
    # DOC: add description 
    """
    Bends atoms 
    w.r.t 0X axis... 
    
    """
    def __new__(self, non_bent_atoms,
                   wedge_angle=None,
                   x_bent_ceter=None,
                   bend_factor=None):
        # FIXME: move(translat) from center mass???
        
        """
        Bends atoms 
        
        Bend_factor definition used:
        
            bend_factor  := (x_max-x_min) / (2 * x_bent_ceter),
        
        which is generalized from bend_factor definition for a tube:
        
            bend_factor  := D / (2*R),
        
        where D is tube diameter and R is the radius of 
        curvature measured from tube axis.  
        
        ...
        """
#        """Either wedge_angle or x_bent_ceter must be given
#        
#        Parameters:
#           
#        @param non_bent_atoms: ase.Atoms 
#            atoms before bent
#        @param wedge_angle: float, must be in range 
#            (WedgeAtoms._min_angle, WedgeAtoms._max_angle)
#            wedge angel 
#            
#        w.r.t 0X axe... 
#        
#            
#        Example 
#        
#        >>> bent_atoms = BendAtoms(atoms, bend_factor=0.95)
#        >>> bent_atoms = BendAtoms(atoms, x_bent_ceter=3.6)
#        >>> bent_atoms = BendAtoms(atoms, wedge_angle=pi/8)
#        >>> bent_atoms = BendAtoms(atoms, wedge_angle=pi/8)
#        w.r.t center mass
#        ... 
#        """
        
        # translate!   
        # ~center = [(x_max - x_min)/2, (y_max - y_min)/2, 0]
        # to be at the center!
        #non_bent_atoms.translate([(x_max - x_min)/2, (y_max - y_min)/2]) 
        
                
        # bent atoms first, then assigne 
        
        # DOC: geometry assumed 
        
        coords = non_bent_atoms.positions.transpose()
        x_min = coords[0].min() 
        x_max = coords[0].max()
        y_min = coords[1].min()
        y_max = coords[1].max()
        
        
        
        # FIXME: geometry/translate!
        # works if only y_min == 0
        if (x_bent_ceter is None) and (wedge_angle is not None):
            # works if only y_min == 0
            x_bent_ceter = y_max / wedge_angle
            bend_factor = (x_max - x_min) / 2 * x_bent_ceter
            
        if (x_bent_ceter is not None) and (wedge_angle is None):
            # works if only y_min == 0
            wedge_angle = y_max / x_bent_ceter
            bend_factor = (x_max - x_min) / 2 * x_bent_ceter
        
        if  ((x_bent_ceter is None) and (wedge_angle is None) 
             and (bend_factor is not None)):
            # works if only y_min == 0
            #bend_factor = (x_max - x_min) / 2 * x_bent_ceter  
            x_bent_ceter  = (x_max - x_min) / 2 * bend_factor 
            wedge_angle = y_max / x_bent_ceter
                
        
        bent_atoms = WedgeAtoms(wedge_angle, non_bent_atoms)
        
        
        # FIXME: use function over list element!
        
        for atom in bent_atoms:
            # FIXME: use numpy!
            atom.position = self._bend_an_atom(self,
                                    atom.position.tolist(),
                                    x_bent_ceter)
              

        #bent_atoms = self._wedge_markup(self, bent_atoms, 
        #                                x_min, x_max, wedge_angle)
        
        return bent_atoms
    
                 
    @staticmethod
    def _bend_an_atom(self, old_coordiantes, x_bent_ceter):
        """
        @param old_coordiantes: list of old coordinates 
        @param x_bent_ceter: 
        """
        
        x_old = old_coordiantes[0]
        y_old = old_coordiantes[1]
        z_old = old_coordiantes[2]
        
        r_new = x_old
        phi_new = y_old / x_bent_ceter
                   
        x_new = r_new * cos(phi_new)
        y_new = r_new * sin(phi_new)
        z_new = z_old
   
        return (x_new, y_new, z_new)
    
    @staticmethod
    def _wedge_markup(self, watoms, x_min, x_max, wedge_angle):
        watoms += Atoms('O', [[x_max * cos(wedge_angle),
                                  x_max * sin(wedge_angle),
                                  0]])
        watoms += Atoms('O', [[x_min * cos(wedge_angle),
                                  x_min * sin(wedge_angle),
                                  0]])
        watoms += Atoms('O', [[x_max * cos(0),
                                  x_max * sin(0),
                                  0]])
        watoms += Atoms('O', [[x_min * cos(0),
                                  x_min * sin(0),
                                  0]])
        return watoms
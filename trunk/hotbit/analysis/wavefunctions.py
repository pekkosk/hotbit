import box
import numpy as nu
import math
from box import mix
acos=math.acos
cos=math.cos
sin=math.sin
atan=math.atan
sqrt=math.sqrt
pi=math.pi


        
def phival(x,y):
    """ Return azimuthal angle for ALL x,y. """
    e=1E-16
    if x>e and y>e:
        return atan(y/x)
    elif x<-e and y>e:
        return atan(y/x) + pi
    elif x<-e and y<-e:
        return atan(y/x) + pi
    elif x>e and y<-e:
        return atan(y/x)+2*pi
    elif abs(x)<=e and abs(y)<=e:
        return 0.0
    elif x>e and abs(y)<=e:
        return 0.0
    elif y>e and abs(x)<=e:
        return pi/2
    elif x<-e and abs(y)<=e:
        return pi
    elif y<-e and abs(x)<=e:
        return 3*pi/2    
    else:
        raise RuntimeError('Strange things in phival')        
        
                
states=['s','px','py','pz','dxy','dyz','dzx','dx2-y2','d3z2-r2']     
def angular(r,wf):
    """ Return angular part of wave function.
    
    The real combinations of Y_lm's
    
    parameters:
    -----------
    r: (not normalized) position vector
    wf: index or symbol for state (look below)
    """    
    R=sqrt(sum(r**2))
    if R<1E-14:
        return 0.0
    theta=acos(r[2]/R)
    phi=phival(r[0],r[1])
    
    if type(wf)!=type(1):
        wf=states.index(wf)
        
    if wf==0:   return 1/sqrt(4*pi)
    elif wf==1: return sqrt(3/(4*pi))*sin(theta)*cos(phi)
    elif wf==2: return sqrt(3/(4*pi))*sin(theta)*sin(phi)
    elif wf==3: return sqrt(3/(4*pi))*cos(theta)
    elif wf==4: return sqrt(15/(4*pi))*sin(theta)**2*cos(phi)*sin(phi)
    elif wf==5: return sqrt(15/(4*pi))*sin(theta)*cos(theta)*sin(phi)
    elif wf==6: return sqrt(15/(4*pi))*sin(theta)*cos(theta)*cos(phi)
    elif wf==7: return 0.5*sqrt(15/(4*pi))*sin(theta)**2*cos(2*phi)
    elif wf==8: return 0.5*sqrt(5/(4*pi))*(3*cos(theta)**2-1)


class WaveFunctions:
    def __init__(self,atoms,mesh=(20,20,20)):
        """ Wave functions into real grid object. 
        
        parameters:
        -----------
        mesh: x,y,and z divisions for the grid
        """
        calc=atoms.get_calculator()
        self.atoms=atoms
        self.calc=calc
        self.el=calc.el
        self.st=calc.st
        cell=self.atoms.get_cell()
        self.L=nu.array([cell[0,0],cell[1,1],cell[2,2]])
        self.N=mesh
        
    def write_vtk(self,i,fname=None):
        """ Write .vtk file of wave function with *index* i. """
        wf=self.st.wf[:,i].copy()
        orbs=self.el.orbitals()
        wfg=nu.zeros(self.N)
        grid=[]
        for i in range(3):
            grid.append( nu.linspace(0,self.L[i],self.N[i]) )
        for orb,c in zip(orbs,wf):
            symb, orbtype, Rnl=orb['symbol'], orb['orbital'], orb['Rnl']
            for i,x in enumerate(grid[0]):
                for j,y in enumerate(grid[1]):
                    for k,z in enumerate(grid[2]):
                        r0=nu.array([x,y,z])
                        r=self.el.vector(orb['atom'],rj=r0)
                        wfg[i,j,k]+=c*Rnl(mix.norm(r))*angular(r,orbtype)
        box.vtk.rectilinear_vtk(grid,wfg,fname)        
                              
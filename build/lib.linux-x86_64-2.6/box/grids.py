import ase
import numpy as np
from box import mix
from box.interpolation import TrilinearInterpolation
from ase.io.cube import read_cube_data
vec=np.array

class GridData:
    def __init__(self, atoms, data):
        """
        
        Paramaters:
        ===========
        atoms:  ase Atoms instance
        data:   grid data (in the cell of atoms, 3-dimensional array) or
                cube file containing the data
        """
        if type(data)==type(''):
            self.data, aux=read_cube_data(data)
        else:
            self.data=data
        self.atoms=atoms
        self.cell = atoms.get_cell()
        assert not self.cell[2, :2].any() and not self.cell[:2, 2].any()
        self.ldos = None
        self.heights = None
        self.axes=vec([self.cell[i,i] for i in range(3)])     
        self.dg=None
        self.ng=self.data.shape
        self.grids=[np.linspace(0,self.axes[i],self.ng[i]) for i in range(3)]
        self.ztol=0.001
        self.dataip=TrilinearInterpolation(self.data,grids=self.grids)
        
        
    def get_grids(self):
        """ Return the grids points in x-,y-, and z-directions. """
        from copy import copy
        return copy(self.grids)
        
       
    def get_data(self):
        """ Return the data as 3-dimensional array. """
        return self.data        
    

    def get_z_averaged_data(self, z):
        """ Get averaged data for given z."""
        sum=0
        for x in self.grids[0]:
            for y in self.grids[1]:
                sum+=self.dataip(vec([x,y,z]))
        c=sum/(self.ng[0]*self.ng[1])
        return c
        
        
    def get_max_over_xy(self, z):
        """ Get maximum data value in the xy-plane at height z."""
        mx=0.
        for x in self.grids[0]:
            for y in self.grids[1]:
                mx=max(mx,self.dataip(vec([x,y,z])))
        return mx
                      
        
    def scan(self,threshold,bottom=0.0):
        """
        Scan the data height profile in constant contour mode with given threshold.
        
        Construct the profile into the same existing x-y grid.
        """
        heights=np.empty((self.ng[0],self.ng[1]))
        nz=len(self.grids[2])
        for i,x in enumerate(self.grids[0]):
            for j,y in enumerate(self.grids[1]):
                # approach from top...
                for k in range(nz-1,-1,-1):
                    if k==0:
                        heights[i,j]=bottom #no crossing of threshold 
                    else:
                        if self.data[i,j,k]<=threshold<=self.data[i,j,k-1]:
                            heights[i,j]=self.find_height(threshold,x,y,self.grids[2][k-1],self.grids[2][k])
                            break
        self.heights=heights   
        return heights      
        
        
    def find_height(self,threshold,x,y,zmin,zmax):
        """
        Let us approach from zmax towards zmin. Return z where the threshold value
        for LDOS is obtained. Search is done recursively. Return NaN
        if threshold not found.
        """
        zmid=0.5*(zmax+zmin)
        if zmax-zmin<self.ztol:
            return 0.5*(zmax+zmin)
        else:
            if self.dataip(vec([x,y,zmid]))>=threshold:
                return self.find_height(threshold,x,y,zmid,zmax)
            else:
                return self.find_height(threshold,x,y,zmin,zmid)
                       
                
    def linescan(self, start, end, threshold, points=200, bottom=0.0):
        """ Make data linescan between start and end points.
        
        Note: only x- and y-coordinates of the start and end points matter.
        
        Parameters:
        -----------
        start: position array
        end: position array
        threshold: constant-data contour line threshold
        points: linear interpolation in output using this many points.
        """
        r1=vec(start)
        r2=vec(end)
        L=np.linalg.linalg.norm(r2-r1)
        nz=len(self.grids[2])
        l=[]
        h=[]
        for i,r in enumerate([r1+i/(points-1.0)*(r2-r1) for i in range(points)]):
            l.append(1.0*i*L/(points-1.0))
            dataz=[ self.dataip(vec([r[0],r[1],self.grids[2][k]])) for k in range(self.ng[2]) ]
            for k in range(nz-1,-1,-1):
                if k==0:
                    h.append(bottom) #no crossing of threshold; use 'bottom' value
                else:
                    if dataz[k]<=threshold<=dataz[k-1]:
                        h.append( self.find_height(threshold,r[0],r[1],self.grids[2][k-1],self.grids[2][k]) )
                        break
        return l,h
                
    def broadened_linescan(self, start, end, width, threshold, points=100, nwidth=10, bottom=0.0, getlines=False):
        """ Make data linescan between start and end points with effective 'width' for data scanner.
        
        Not: only x- and y-coordinates of the start and end points matter.
        Line scan between start and end, but with such a way that you perform
        'nwidth' parallel line scans from width/2 distance from both sides of the original line.
        The result with given distance from the start is the maximum over all parallel values.
        
        Parameters:
        -----------
        start: position array
        end: position array
        width: width (in Angstrom) for data scanning 
        threshold: constant-data contour line threshold
        nwidth: how many parallel line scans are performed
        points: linear interpolation in output using this many points.
        getlines: if True, return ase.Atoms object where atoms are 'X's with 
                  starting and end-points of the line scans
                  
        Return: distance array and heights on that distance array
                (if getlines, return also the ase.Atoms object)
        """  
        start, end = vec(start), vec(end)
        n=np.cross( end-start,vec([0,0,1]) )
        n/=np.linalg.norm(n)             
        start0, end0 = start - (width/2)*n, end - (width/2)*n
        atoms=ase.Atoms()
        
        # perform nwidth parallel scans        
        hlist=[]
        for shift in np.linspace(0,width,nwidth):
            l,h = self.linescan(start0+shift*n, end0+shift*n, threshold, points=points, bottom=bottom)
            atoms.append( ase.Atom('X',start0+shift*n) )
            atoms.append( ase.Atom('X',end0+shift*n) )
            hlist.append(h)                    
        hlist=np.array(hlist)                    
        
        # take maxima from these scans
        h=[]
        for i in range(len(l)):
            h.append( max(hlist[:,i]) )
        if getlines:            
            return l,h,atoms            
        else:
            return l,h
            
            
            
       
    

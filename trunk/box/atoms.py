# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

import numpy as npy
from box import mix
from box.data import data

from ase import units
from ase import Atom as ase_Atom
from ase import Atoms as ase_Atoms
from copy import copy
vec=npy.array


class Atoms(ase_Atoms):
    """A class for cluster structures
    to enable simplified manipulation"""

    def __init__(self, *args, **kwargs):


        if kwargs.get('filename') is not None:
            filename = kwargs.pop('filename')
            ase_Atoms.__init__(self, *args, **kwargs)
            #ase_Atoms.read(self,filename, kwargs.get('filetype'))
        #elif isinstance(args[0],ase_Atoms):
            #self=args[0]
            #print 'dsdsdsdsdgsdgsdgsdg'
            #assert False
        else:
            ase_Atoms.__init__(self, *args, **kwargs)
            
        self.data={}
        self.set_data()
        self.analyzed=[]
        
    def set_calculator(self,calc=None):
        ase_Atoms.set_calculator(self,calc)
        if calc!=None:
            self._update()
                    
    def set_positions(self,r):
        """ If calc is set, system should be reasonable for updating. """
        ase_Atoms.set_positions(self,r)
        if self.calc!=None:
            self._update()
        
    def set_cell(self,cell,**kwargs):
        """ If calc is set, system should be reasonable for updating. """
        ase_Atoms.set_cell(self,cell,**kwargs)
        if self.calc!=None:
            self._update()        
        
    def set_data(self):
        """ Read element data in. """
        symbols=self.get_chemical_symbols()
        keys=data[symbols[0]].keys()
        for key in keys:
            self.data[key]=[data[symbol][key] for symbol in symbols]
        
    def get_name(self):
        """ Get the name of the system, e.g. 'benzene'. Default are symbols, 'H2O' """
        if 'atoms_name' in self.data:
            return self.data['atoms_name']
        else:
            symbols=self.get_chemical_symbols()
            dict={}
            for symbol in symbols:
                n=symbols.count(symbol)
                if n not in dict:
                    dict.update({symbol:n})
            name=''
            for symbol in dict:
                if dict[symbol]<2:
                    name+=symbol
                else:
                    name+=symbol+str(dict[symbol])
            self.data['atoms_name']=name
            return name
                
    def get_N(self):
        return len(self.get_positions())
            
    def set_name(self,name):
        """ Give the system a name, e.g. 'benzene' """
        self.name=name
        
    def _update(self):
        """ After changing cell or atom positions, make checks etc. """
        self._check_cell_orthogonality()
        self._reduce_atoms_into_cell()
        
    def _reduce_atoms_into_cell(self):
        """
        Collect all atoms into the 'first' unit cell (0<=x<L_x,...)
        """
        L=self.get_cell_axes()
        r=self.get_positions()
        if npy.any([x<1.1 for x in L]):
            raise ValueError('Cell (%f,%f,%f) too small and not set yet.' %(L[0],L[1],L[2]))
        for i in range(len(r)):
            for a in range(3):
                if self.pbc[a]:
                    r[i,a]=npy.mod(r[i,a],L[a])
                elif r[i,a]<0 or r[i,a]>L[a]:
                    raise ValueError('Atom %i at %f outside box in non-pbc direction %i' %(i,r[i,a],a) )
        ase_Atoms.set_positions(self,r)        
    
    def _check_cell_orthogonality(self):
        """
        Non-orthogonality of unit cell is not implemented yet.
        """
        for a in range(3):
            aux=range(3)
            aux.pop(a)
            if npy.any(self.cell[a,aux]>1E-9): 
                raise NotImplementedError('Unit cell not orthogonal.')               
        return True
                
    def get_cell_axes(self):
        """
        Lengths of the unit cell axes (currently for orthonormal cell).
        """
        return self.cell[0,0],self.cell[1,1],self.cell[2,2]
        
    def vector(self,ri,rj=None,mic=True):
        """ 
        Return the vector ri or rj-ri, normally within the minimum image convention (mic).
        
        parameters:
        -----------
        ri: initial point (atom index or position)
        rj: final point (atom index or position)
        """
        if not mic:
            raise NotImplementedError('Why not use mic?')
        R=self.positions
        if type(ri)==type(1):
            Ri=R[ri].copy()
        else:
            Ri=ri.copy()
        if rj==None:
            Rj=None
        elif type(rj)==type(1):
            Rj=R[rj].copy()
        else:
            Rj=rj                                    
            
        if Rj==None:
            rij=Ri
        else:
            rij=Rj-Ri
            
        L=self.get_cell_axes()
        for a in range(3):
            if self.pbc[a]:
                rij[a]=npy.mod(rij[a]+0.5*L[a],L[a]) - 0.5*L[a]
        return vec(rij)
            
    def distance(self,ri,rj=None,mic=True):
        """
        Return the length of |ri| or |rj-ri|, normally within the minimum image convention.
        """
        if rj==None:
            return mix.norm(self.vector(ri,mic))
        else:
            return mix.norm(self.vector(ri,rj,mic))
    
    def displace_atoms_randomly(self,dR=0.05):
        """
        Displace atoms with random positions up to given radii 
        with Gaussian profile (sigma=dR/2).
        Same for all atoms or list of radii for all atoms.
        """
        from numpy.random import normal
        sin=npy.sin
        cos=npy.cos
        R=self.get_positions()
        if type(dR)!=type([]):
            dR=[dR]*len(R)
        for i in range(len(R)):
            phi,theta=mix.random_direction()
            disp=min(normal(scale=dR[i]/2),dR[i])
            R[i]=R[i]+disp*vec([sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)])
        self.set_positions(R)
                  
    def scale_positions(self,x):    
        """ Scale the whole system by x; also the unit cell. """
        self.set_cell(self.get_cell()*x,scale_atoms=True)
        self._reduce_atoms_into_cell()    
        
        
    def copy(self):
        """ Override ase atoms copy (not to get ase.Atoms type) """        
        from copy import copy
        return copy(self)
        #atoms = Atoms(cell=self.cell, pbc=self.pbc)
        #atoms.arrays = {}
        #for name, a in self.arrays.items():
            #atoms.arrays[name] = a.copy()
        #return atoms
        
    def spatial_spans(self):
        r=self.get_positions()        
        return vec([max(r[:,a])-min(r[:,a]) for a in range(3)])    
    
    def get_kinetic_temperature(self,fixed=0,internal=False):
        """ Use equipartition to get temp (in K) without 'fixed' number degr. of freedom. """
        if not internal:
            return self.get_kinetic_energy()*2 / ((3*self.get_N()-fixed)*units.kB)
        elif internal:
            assert fixed==0
            raise NotImplementedError('todo: get kinetic energy with removed cm translation and rotation')
    
    def get_total_energy(self):
        """ Return potential+kinetic energy. """
        return self.get_kinetic_energy() + self.get_potential_energy()
        
    #def __iter__(self):
        #""" Initialize iterator. """
        #self.itr=-1
        #return self
            
    #def next(self):
        #""" Iterator over atoms. """
        #self.itr+=1
        #if self.itr==self.N: 
            #raise StopIteration
        #return self.atoms[self.itr]        
            
    def set_cell(self,cell,scale_atoms=False):
        """ Override Atoms.set_cell with default atom fixing. """
        ase_Atoms.set_cell(self,cell,scale_atoms=scale_atoms)
        
    def angular_momentum(self):
        assert NotImplementedError()
        
    def moments_of_inertia(self,rcm=True):
        """ Return the moments and principal axes of inertia (wrt. center of mass if rcm is True) """        
        r=self.get_positions()
        mass=self.get_masses()
        N=len(r)
        if rcm:
            cm=self.get_center_of_mass()
            r-=cm            
            
        # make inertia tensor and solve its eigenvalues (principal moments
        # of inertia) and eigenvectors (principal axis of inertia)
        tensor=npy.zeros((3,3))    
        diag=sum([mass[k]*sum(r[k,:]**2) for k in range(N)])
        for i in range(3):
            tensor[i,i]=diag
            for j in range(3):
                tensor[i,j]=tensor[i,j]-sum([mass[k]*r[k,i]*r[k,j] for k in range(N)])
                    
        eigs,vectors=npy.linalg.eigh(tensor)
        return eigs,vectors
        
    def center_of_mass_to_origin(self):
        assert NotImplementedError()
        
    def remove_rotations(self):
        assert NotImplementedError()
        
    def atom_distribution_wrt_point(self,point=None):
        assert NotImplementedError()
        #of point==None, point=center_of_mass
                
    def spatial_extrema(self):
        assert NotImplementedError()
        #return minx,maxx,miny,maxy,..
    
    def radial_distribution_function(self):
        assert NotImplementedError()
        #within minimum image convention (max of function= minimum box length/2)
    
    def coordination_numbers(self):
        assert NotImplementedError()
        
    def mean_bond_length(self):
        """ Return mean bond length using pair distribution function. """
        pdf = self.pair_distribution_function()
        max = False
        eps = 1e-6
        # find the first max and following min of pair distr. fct
        if self.get_N()==2:
            return self.distance(0,1)
        else:
            for i in xrange(len(pdf)-1):
                if pdf[i,1]>pdf[i+1,1]: max=True
                if max==True and (pdf[i,1]<pdf[i+1,1] or pdf[i,1]<eps or i==len(pdf)-3):
                    rcut = pdf[i,0]
                    rij  = self.pair_distribution_list()
                    return sum(npy.where(rij<=rcut,rij,0.0)) / sum(rij<=rcut)
    
    def number_of_bonds(self):
        """ Return the number of bonds (intended for homonuclear molecule). """
        pdf = self.pair_distribution_function()
        max = False
        eps = 1e-6
        # find the first max and following min of pair
        # distr. function
        for i in xrange(len(pdf)-1):
            if pdf[i,1]>pdf[i+1,1]: max=True
            if max==True and (pdf[i,1]<pdf[i+1,1] or pdf[i,1]<eps):
                break
        rcut = pdf[min(i,len(pdf)),0]
        rij  = self.pair_distribution_list()
        return sum( rij<=rcut )
            
    def pair_distribution_function(self,rmin=0.0,rmax=10,sigma=0.7):
        """ Return Gaussian-broadened pair distribution function (or RDF). """
        if rmax==None:
            rmax=max(min(self.get_cell_axes())/2,5)
        rij=self.pair_distribution_list()
        rij.sort()
        grd=npy.linspace(rmin,rmax,200)
        pdf=vec([grd,npy.zeros(len(grd))]).transpose()
        for x in rij:
            for i in xrange(len(grd)):
                pdf[i,1]+=mix.gauss_fct(grd[i],mean=x,sigma=sigma)
        return pdf
        
    def pair_distribution_list(self):
        """ Return the array r_ij=|r_i-r_j| for all pairs. """
        r=self.get_positions()
        N=len(r)
        rij=[]
        for i in xrange(N):
            for j in xrange(i+1,N):
                rij.append( self.distance(i,j) )
        return vec(rij)        
            
    def list_properties(self):
        """ Return all the data properties """
        return self.data.keys()
    
    def get(self,key):
        if key in self.data:
            return self.data[key]
        else:
            return None       
            
    def set(self,key,value):
        ind=1
        try:
            ind=len(value.shape)
        except:
            ind=1
        if ind==2:
            self.data[key]=[value[i,:] for i in range(len(value[:,0]))]
        else:
            self.data[key]=value
        
    def construct_bonds(self):
        """
        Make the bond list for the molecule.
        
        Use estimates for bond lengths from van der Waals radii.
        Make bond if R<(R_cov,1+R_cov,2)*1.4
        """
        nb=0
        bonds=[]
        r=self.get_positions()
        rc=self.get('R_cov')
        for i in range(len(r)):
            for j in range(i+1,len(r)):
                R=self.distance(i,j)
                if R<(rc[i]+rc[j])*1.4: 
                    bonds.append({'i':i,'j':j,'length':R})
                    nb+=1
        self.set('nbonds',nb)
        self.set('bonds',bonds)
        
        
    def write_vtk(self,file=None):
        """ vtk output of atoms (all scalar and vector properties) """
        if file==None:
            file=self.get_name()+'.vtk'
        N=len(self)
        f=open(file,'w')
        f.write('# vtk DataFile Version 2.0 \nAtoms %s\n' %self.get_name())
        f.write('ASCII \nDATASET UNSTRUCTURED_GRID\n')
        f.write('POINTS %i double \n ' %N)
        fmt='%20.14f' #output format for floats
        
        # Point data (atom coordinates) and cell data (bonds)
        if self.get('nbonds')==None: 
            self.construct_bonds()
        nb=self.get('nbonds')
        bonds=self.get('bonds')
        for r in self.get_positions():
            f.write('%s\n' %mix.a2s(r,fmt=fmt))
        f.write('CELLS %i %i\n' %(nb,3*nb))
        for bond in bonds:
            f.write( '2 %i %i\n' %(bond['i'],bond['j']) )
        f.write('CELL_TYPES %i\n' %nb)
        for bond in bonds:
            f.write('3\n')    
            
        # First the data related to atoms
        f.write('POINT_DATA %i\n' %N)
        self.set('velocities',self.get_velocities())        
        for property in self.list_properties():
            properties=self.get(property)
            try:
                tp=type(properties[0])
            except:
                continue            
            if tp==type(vec([])) or tp==type([]):
                if not isinstance(properties[0],(int,float)): continue
                print>>f, 'VECTORS %s double\n' %property
                for value in properties:
                    print>>f, mix.a2s(value,fmt=fmt)
            else:
                try:
                    x=float(properties[0])
                except:
                    continue                  
                print>>f, 'SCALARS %s double 1\nLOOKUP_TABLE default' %property
                for value in properties:
                    print>>f, '%12.6f' %(value*1.0)
            
            
        # Then data related to bonds
        print>>f, 'CELL_DATA %i' %nb
        print>>f, 'SCALARS bond_length double 1\nLOOKUP_TABLE default'
        for bond in bonds:
            f.write( '%f\n' %bond['length'] )
        f.close()            

        
        
        
    #def __add__(self, other):
        #""" 
        #Add (one) atom or another molecule. 
        #Note: completely new molecule constructed.
        #"""
        #new=copy(self.atoms)
        #if type(other)==type(Molecule()):
            #for atom in other:
                #new.append(atom)
        #else:
            #new.append(other)
        #return Molecule(new)
        
     
     

    
    #def get_charge(self):  
    #def copy to periodic directions (gpaw->hotbit calculations...)  
    #def get(self,property):
    
    
    
    #def sp_occupations(self,excess_el=0):      
         #energy_as_separated(self,excess_el=0):
        #return IE,EA
                                     
                

                                                         
    #def read_bonds(self,file):
        #"""
        #Read bond information from given file.
        
        #* file -- file name of file object where the bonding info is read (for
          #file object the "next" data set is read.)
        #"""
        #self.N_bonds=0
        #f,opened=mix.file_safeopen(file)
        #lines=mix.read(file,fmt='string')
        #self.bonds=[]
        #print file
        #for line in lines:
            #self.N_bonds +=1
            #i,j,B,r = line.split()
            #self.bonds.append( [int(i)-1,int(j)-1,float(B)] )                                                          
            
    
    
    

      
                

if __name__=='__main__':
    print 'atoms test run'
    atoms=Atoms(symbols='H2O',positions=[(1,0,0),(0,1,0),(0,0,0)])
    atoms.set_cell([50,50,50])
    atoms.set_pbc(True)
    print 'H2O:\n',atoms.get_positions()
    atoms.reduce_atoms_into_cell()
    print 'H2O:\n',atoms.get_positions()
    print atoms.vector(vec([0,0,0]),vec([0,4.9,0]))
    
    print 'mean bond length',atoms.mean_bond_length()
    print 'number of bonds',atoms.number_of_bonds()
    #print 'pair list',atoms.pair_distribution_list()
    print atoms.get_name()
    atoms.write_vtk('koe.vtk')
    
    #print 'pdf',atoms.pair_distribution_function()
    #pdf=atoms.pair_distribution_function()
    #import pylab as p
    #p.plot(pdf[:,0],pdf[:,1])
    #p.show()
    ##many=Atoms(symbols='H500',positions=[(10,10,10) for x in range(500)],cell=(20,20,20),pbc=True )
    #many.displace_atoms_randomly(dR=7.0)
    #from ase import write
    #write('shake.xyz',many)
    

    

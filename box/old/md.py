#!/usr/bin/env python
"""
    Contain classes related to molecular calculations. 
    
    Atomic units used throughout unless otherwise stated.
    
    Author P. Koskinen 15.9 2006
    
"""
import os,sys
import box.mix as mix
from numpy import *

class Element:
    """
    Class for elements. 
    
    Contains element-specific information such
    as mass, name, ionization potentials etc.
    """
    def __init__(self,element,fil=None):
        """
        Initialize the element object.
        
        At least symbol given, and if file (=elements.dat used for HOTBIT)
        is present, read more element info from there.
        
        Parameters
            
            * element -- element symbol
            
            * fil -- the path of elements.dat -file
            
        """
        self.element=element.strip()
        if fil:
            self.read_element_info(fil)
               
    def __call__(self):
        """
        Print the element data (__dict__)
        """
        print self.__dict__
    
    def read_element_info(self,file):
        """
        Read data for the element from elements.dat-file.
        
        * file -- elements.dat -file
        """
        f=open(file)
        lines=f.readlines()
        # read until the beginning of the element
        for i in range( len(lines) ):
            data=lines[i].split('=')
            if len(data)<=1: continue # not a data line
            if data[0].strip()=='element':
                hlp=data[1].split()
                if hlp[0].strip()==self.element: break
        
        # read the data for the element
        start=i+1
        if start>=len(lines): 
            mix.error_exit('No element '+self.element)
        for i in range( start,len(lines) ):
            data=lines[i].split('=')
            if len(data)<=1: continue # not a data line
            key=data[0].strip()
            if key=='element': break # don't read new element
            val=data[1].strip()
            if key=='common' : self.common=val
            if key=='Z'      : self.Z=int(val)
            if key=='m'      : self.m=float(val)
            if key=='q0'     : self.q0=float(val)
            if key=='no'     : 
                self.no=int(val)
                self.oe=[float(x) for x in lines[i+1].split()]
            if key=='U'      : self.U=float(val)
            if key=='FWHM'   : self.FWHM=float(val)
            if key=='vib'    : self.vib=float(val)      
            if key=='dnn'    : self.dnn=float(val)      
            if key=='IE'     : self.IE=float(val)
            if key=='EA'     : self.EA=float(val)
            if key=='R_vdw'  : self.R_vdw=float(val)
            if key=='R_cov'  : self.R_cov=float(val)
        f.close()
    
    def sp_occupations(self,excess_el=0):
        """
        Return the occupations for single particle states with given N_el.
        
        * excess_el -- number of excess electrons on element (wrt. neutral)
        """
        q0 = self.q0 + excess_el
        N_full, rem = q0//2., q0%2.
        f=zeros(self.no)
        f[0:N_full]=2.0
        f[N_full]  =rem
        return f 
       
    def energy_as_separated(self,excess_el=0):
        """
        Return the energy of isolated atom .
        
        Return sum_i occ_i e_i.
        
        * excess_el -- number of excess electrons on element (wrt. neutral)
        """
        en    = sort( self.oe )
        f     = self.sp_occupations(excess_el)
        ebs   = sum( [e*f for e,f in zip(en,f)] )
        ecoul = 0.5*self.U*excess_el**2
        return ebs+ecoul
    
    def energies(self,ret):
        """
        Return some characterisics of the electronic structure.
        
        * ret -- the returned energy
        
            * IE -- ionization energy
            * EA -- electron affinity
        """
        E=self.energy_as_separated
        if ret=='IE':
            return ( E(-1)-E(0) )
        elif ret=='EA':
            return ( E(0)-E(1) )
   
class Atom:
    """
    Class for atoms. 
    
    Atom consists the element information, using the class
    'Element', but contains also position, velocity etc. variable information.
    """
    def __init__(self,elem,r=array([0.,0.,0.]),v=array([0.,0.,0.]),\
                            f=array([0.,0.,0.])):
        """
        Initialize atom information.
        
        * elem -- the element of atom
        * r -- position-vector
        * v -- velocity-vector
        * f -- force-vector
        """
        self.el=elem
        self.r=r
        self.v=v
        self.f=f
        self.flag=1
        self.dq=0.0 # atom additional electron population
        self.vectors={} # dictionary of general vector data for atoms
        
      
        
        
class Molecule:
    """
    Class for molecules. 
    
    Consists of many atoms (using class 'Atom'), and
    hosts also other additional information such as number of atoms,
    electrons in the whole molecule, binding info.
    """
    def __init__(self,file=None,format='xyz',efile=None):
        """
        Initialize molecule.
        
        * file -- if present, read molecule from this file
        * format -- format of the given file ('xyz','I_info', or 'dat')
        * efile -- the path for elements.dat -file used for element info
        """
        from os import environ
        self.atoms=[]
        self.vector_list=[] # list of names of vectors defined for all atoms
        self.N=0        # number of atoms in molecule
        self.vel=0      # number of valence electrons
        self.N_bonds=0  # number of bonds
        if efile:   
            self.efile=efile            
        else:       
            self.efile=environ.get('BOX_ELEMENTS')
            if not self.efile: mix.error_exit('elements.dat file not specified')
        if file:
            if format=='xyz'      : self.input_xyz(file)
            if format=='I_info'   : self.input_I_info(file)
            if format=='dat'      : self.input_atoms_dat(file)
               
            
    def add_molecule(self,mol2):
        """ 
        Adds the atoms of another molecule to the present one.
        
        * mol2 -- the molecule the atoms of which will be added
        """
        for atom in mol2.atoms:
            self.add(atom)

    def add(self,atom):
        """
        Add atom into the molecule.
        
        * atom -- the atom to be added 
        """
        self.N+=1
        self.atoms.append(atom)
        
    def input_atoms_dat(self,file):
        """
        Read molecule from atoms.dat-file.
        
        * file -- the given atoms.dat -file
        """
        fnd=mix.find_value
        f=open(file,'r')
        N=int(fnd(f,'nat'))
        self.extra_electrons=float(fnd(f,'extra_electrons',default=0.0))
        atoms=fnd(f,'atoms',fmt='strings')
        for line in atoms:
            el,x,y,z = line.split()
            r=array([float(x),float(y),float(z)])
            elem=Element(el,self.efile)
            atom=Atom(elem,r)
            self.add( atom )
            
        if N!=self.N: error_exit('Error in atoms.dat file')
            
        flags=fnd(f,'flags',fmt='strings',default='def')
        if flags!='def': 
            for i in range(self.N): self.atoms[i].act=int(flags[i])
        else:            
            for i in range(self.N): self.atoms[i].act=1
        vels=fnd(f,'velocities',fmt='strings',default='def')
        if vels!='def': 
            for i in range(self.N):
                self.atoms[i].v=array([float(x) for x in vels[i].split()]) 
        
    def input_I_info(self,file):
        """
        Reads molecule from I_info -file, used for Cmdft program.
        
        * file -- the I_info file for input
        """
        fi=open(file)
        flines=fi.readlines()
        N=0
        for line in flines:
            if line.find("\"")>=0: N+=1
        for i in range(N):
            line=flines[i+1]
            symb=line[1:3]
            rest=line[4:]
            m,vel,x,y,z=rest.split()
            vx,vy,vz=flines[N+(2*i+1)].split()
            fx,fy,fz=flines[N+(2*i+2)].split()
            r=array([float(x ),float(y) ,float(z)])
            v=array([float(vx),float(vy),float(vz)])
            f=array([float(fx),float(fy),float(fz)])
            elem=Element(symb,self.efile)
            elem.m=float(m)
            elem.vel=float(vel)
            self.add( Atom(elem,r,v,f) )
        fi.close()
            
            
    def input_xyz(self,file):
        """
        Read molecule from xyz-file. 
        
        * file -- the given input xyz-file. It can be either a file name
        or a file object. If file is a file object the next "frame" is read
        without rewinding the file.
        """
        if isinstance(file,str):
            # input was a file name: read it
            if not os.path.exists(file):
                mix.error_exit('File '+file+' does not exist. Exit now.')
            f=open(file)
            flines=f.readlines()
            f.close()
        else:
            # input was file object: read the next 'frame'
            flines=[]
            next=file.readline()
            if len(next)==0: return
            flines.append( next )
            flines.append( file.readline() )
            N=int(flines[0])
            for i in range(N):
                flines.append( file.readline() )
        for line in flines[2:len(flines)]:
            # read the atoms symbols and coordinates
            try:
                lst = line.strip().split()
                symb,x,y,z=(lst[0],lst[1],lst[2],lst[3])
                #if len(lst)==4: symb,x,y,z=lst
                #else: symb,x,y,z,rest=lst
            except:
                break #probably end of file
            r=array([float(x),float(y),float(z)])/0.529177 #input in A
            elem=Element(symb,self.efile)
            atom=Atom(elem,r)
            self.add(atom)

            
    def output_xyz(self,file):
        """
        Write the molecule into a xyz-file. 
        
        * file -- the output file name or file object. If it is an file
        object, the molecule is simpy appended as next "frame".
        """
        try:
            flag=0
            f=open(file,'w')
        except:
            flag=1
            f=file
        f.write(str(self.N)+'\n')
        f.write('xyz generated by write_xyz [md.py]\n')
        for atom in self.atoms:
            r = atom.r*0.529177
            f.write( '%s %f %f %f\n' %(atom.el.element,r[0],r[1],r[2])  )
        if flag==0: f.close()
        
    def output_I_info(self,file):
        """
        Write the molecule into a I_info-file for Cmdft.
         
         * file -- the output file name.
        """
        f=open(file,'w')
        f.write('0\n')
        for a in self.atoms:
            # print coordinates
            #f.write( "\""+symb+"\" "+vel+" "+m+" "+r[1:-1]+'\n'  )
            f.write( "\"%s\" %6.2f %9.3f %.12g %.12g %.12g\n" \
                    %(a.el.element,a.el.vel,a.el.m,a.r[0],a.r[1],a.r[2]) )
        for a in self.atoms:
            # print forces and velocities
            f.write( "%.12g %.12g %.12g\n" %(a.v[0],a.v[1],a.v[2]) )
            f.write( "%.12g %.12g %.12g\n" %(a.f[0],a.f[1],a.f[2]) )
        f.close()
        
       
    def output_atoms_dat_old(self,file):
        """
        Writes the molecule into an atoms.dat file.
        
        This is for the older (<5.3 2007) version of atoms.dat.
        
        * file -- output file name
        """
        f=open(file,'w')
        f.write("<----------------Total number of atoms\n")
        f.write( str(self.N)+"\n" )
        f.write("<----------------Number of occupied orbitals\n")
        f.write(str(self.vel/2.)+'\n')
        f.write("<----------------Element,mass,coord,act,langevin\n")
        for atom in self.atoms:
            m = str(atom.el.m)
            r = str(atom.r)
            f.write( atom.el.element+" "+m+" "+r[1:-1]+" 1 0 0\n")
        f.write("<----------------velocities\n")
        for atom in self.atoms:
            f.write( str( atom.v )[1:-1]+"\n")
        f.write("<----------------forces\n")
        for atom in self.atoms:
            f.write( str( atom.f )[1:-1]+"\n")
        f.close()
        
    def output_atoms_dat(self,file,extra_electrons=0):
        """
        Write the molecule into atoms.dat-file.
        
        * file -- output file name
        """
        f=open(file,'w')
        f.write("nat=%i\n" %self.N)
        if extra_electrons!=0: f.write('extra_electrons=%i\n' %extra_electrons)
        f.write("atoms=\n")
        for atom in self.atoms:
            x,y,z=atom.r[:]
            el=atom.el.element
            f.write( '%s %f %f %f\n' %(el,x,y,z))
             
        f.close()
        #These are rarely needed, so don't add them (yet)
        #f.write('velocities=\n')
        #for atom in self.atoms: f.write( str( atom.v )[1:-1]+"\n")
        #f.write('forces=\n')
        #for atom in self.atoms: f.write( str( atom.f )[1:-1]+"\n")
            
    def r_cm(self):
        """
        Return the center of mass vector of the molecule.
        """
        r=array([0,0,0])
        m_tot=0
        for atom in self.atoms:
            m=atom.element.m
            m_tot+=m
            r=r+atom.r
        return r/m_tot
        
    def translate(self,v):
        """
        Translate the molecule.
        
        * v -- translate molecule by v
        """
        for atom in self.atoms:
            atom.r=atom.r+v
        
    def scale_r(self,x):
        """
        Scale all the coordinates.
        
        * x -- scaling factor.
        """
        for atom in self.atoms:
            atom.r=atom.r*x
            
    def __call__(self):
        """
        Print some molecule data.
        
        Prints locations.
        """
        print "Molecule:", self.N
        for atom in self.atoms:
            print "%10.4f %10.4f %10.4f" %tuple(atom.r[:])
        
    def r2array(self):
        """
        Return the atom locations as an array for faster 
        computations.
        """
        r=[]
        for atom in self.atoms:
            r.append( atom.r )
        return array(r)

    def pair_distr_list(self):
        """
        Return the array r_ij=|r_i-r_j| for all pairs a 1-D array.
        """
        r=self.r2array()
        rij=[]
        for i in xrange(self.N):
            for j in xrange(i+1,self.N):
                rij.append( sqrt(sum((r[i]-r[j])**2)) )
        return array(rij)    
            
    def pair_distr_function(self,rmin,rmax,sigma=0.7):
        """
        Return the pair distribution function.
         
        * rmin -- the minimum of the function
        * rmax -- the maximum of the function
        * sigma -- Gaussian used in the broadening, this is its sigma
        """
        rij=self.pair_distr_list()
        rij=sort(rij)
        grd=mix.grid(rmin,rmax,200)
        g  =array([grd,zeros(len(grd))]).transpose()
        for x in rij:
            for i in xrange(len(grd)):
                g[i,1]+=mix.gauss_fct(grd[i],mean=x,sigma=sigma)
        return g
        
    def av_bond_length(self):
        """
        Return the average bond length using estimates
        from pair distribution function.
        """
        g   = self.pair_distr_function(0.0,15.0,sigma=0.7)
        max = False
        eps = 1e-6
        # find the first max and following min of pair
        # distr. function
        if self.N==2:
            r=self.r2array()
            return mix.norm(r[:,0]-r[:,1])
        else:    
            for i in xrange(len(g)-1):
                if g[i,1]>g[i+1,1]: max=True
                if max==True and (g[i,1]<g[i+1,1] or g[i,1]<eps or i==len(g)-3):
                    rcut = g[i,0]
                    rij  = self.pair_distr_list()
                    return sum( where(rij<=rcut,rij,0.0) )/sum( rij<=rcut )
        
    def nr_bonds(self):
        """
        Return the number of bonds (for homonuclear molecule).
        """ 
        g   = self.pair_distr_function(0.0,10.0,sigma=0.7)
        max = False
        eps = 1e-6
        # find the first max and following min of pair
        # distr. function
        for i in xrange(len(g)-1):
            if g[i,1]>g[i+1,1]: max=True
            if max==True and (g[i,1]<g[i+1,1] or g[i,1]<eps):
                rcut = g[i,0]
                rij  = self.pair_distr_list()
                return sum( rij<=rcut )
     
        
    def move_atom(self,atom,dr):
        """
        Translate atom by the vector dr.
        
        * atom -- the atom index to be translated. First atom=0.
        * dr -- the translation vector
        """
        self.atoms[atom].r+=dr
         
    def energy_atoms_separated(self):
        """
        Return the total energy of separate atoms.
         
        Sum of sp energies.
        """
        sm=0.
        for atom in self.atoms:
            sm+=atom.el.energy_as_separated()
        return sm
        
    def construct_bonds(self):
        """
        Make the bonding list for the molecule.
        
        Use estimates for bond lengths from van der Waals radii.
        Make bond if R<(R_cov,1+R_cov,2)*1.4
        """
        nb=0
        self.bonds=[]
        for i in range(self.N):
            for j in range(i+1,self.N):
                ai = self.atoms[i]
                aj = self.atoms[j]
                r1 = ai.el.R_cov
                r2 = aj.el.R_cov
                R  = mix.norm( ai.r-aj.r )                
                if R<(r1+r2)*1.4: 
                    self.bonds.append( [i,j,R] )
                    nb+=1
        self.N_bonds = nb
        
    def read_bonds(self,file):
        """
        Read bond information from given file.
        
        * file -- file name of file object where the bonding info is read (for
          file object the "next" data set is read.)
        """
        self.N_bonds=0
        f,opened=mix.file_safeopen(file)
        lines=mix.read(file,fmt='string')
        self.bonds=[]
        print file
        for line in lines:
            self.N_bonds +=1
            i,j,B,r = line.split()
            self.bonds.append( [int(i)-1,int(j)-1,float(B)] )        
    
    def vtk_output(self,fn):
        """
        Make a vtk-file of the current molecule.
        
        Output of coordinates, charges, velocities, forces, bonds, etc.
        
        * fn -- the output file name (e.g. 'molecule.vtk')
        """
        f=open(fn,'w')
        f.write('# vtk DataFile Version 2.0 \nNOTB simulation\n')
        f.write('ASCII \nDATASET UNSTRUCTURED_GRID\n')
        f.write('POINTS %i double \n ' %self.N)
        # 
        # Point data (atom coordinates) and cell data
        # (bonds, connections between points)
        #
        if self.N_bonds==0: self.construct_bonds()
        for atom in self.atoms:
            r=atom.r
            f.write('%f %f %f\n' %(r[0],r[1],r[2]) )
        f.write( 'CELLS %i %i\n' %(self.N_bonds,3*self.N_bonds) )
        for bond in self.bonds:
            f.write( '2 %i %i\n' %(bond[0],bond[1]) )
        f.write( 'CELL_TYPES %i\n' %self.N_bonds )
        for bond in self.bonds:
            f.write( '3\n' )    
           
        #
        # First the data related to atoms
        #    
        f.write('POINT_DATA %i\n' %self.N)
        
        # vdW radii
        f.write('SCALARS R_vdw double 1\nLOOKUP_TABLE default\n')       
        for atom in self.atoms:
            f.write('%f\n' %atom.el.R_vdw)
         
        # element numbers
        f.write('SCALARS Z double 1\n') 
        f.write('LOOKUP_TABLE default\n')       
        for atom in self.atoms:
            f.write('%i\n' %atom.el.Z)
            
        #f.write('COLOR_SCALARS Z 1\n')
        #for atom in self.atoms:
            #c=float(atom.el.Z)/Zm
            #f.write('%f\n' %c)
                        
        #f.write('LOOKUP_TABLE atom_color_table %i\n' %self.N)
        #colors={}
        #colors['H']='1.0 1.0 1.0 1.0'
        #colors['C']='0.0 0.52 0.22 1.0'
        #for atom in self.atoms:
            #f.write( colors[atom.el.element]+'\n' )
    
        #f.write('LOOKUP_TABLE atom_colors 2\n' )
        ##colors={}
        ##colors['H']='0.0 1.0 1.0'
        ##colors['C']='1.0 0.0 1.0'
        ##f.write( colors['H']+'\n' )
        ##f.write( colors['C']+'\n' )
        #f.write('1.0 1.0 1.0 1.0\n')
        #f.write('0.0 0.0 0.0 1.0\n')
        
        #flags
        f.write('SCALARS flag double 1\nLOOKUP_TABLE default\n')       
        for atom in self.atoms:
            f.write('%f\n' %atom.flag)
       
        #charges
        f.write('SCALARS dq double 1\nLOOKUP_TABLE default\n')       
        for atom in self.atoms:
            f.write('%f\n' %atom.dq)
        
        #velocities
        f.write('VECTORS velocities double\n')       
        for atom in self.atoms:
            v=atom.v
            f.write('%f %f %f\n' %(v[0],v[1],v[2]))
                            
        #velocities
        f.write('VECTORS forces double\n')       
        for atom in self.atoms:
            F=atom.f
            f.write('%f %f %f\n' %(F[0],F[1],F[2]))  
              
        #other vector data
        for vec in self.vector_list:
            f.write('VECTORS %s double\n' %vec)       
            for atom in self.atoms:
                v=atom.vectors[vec]
                f.write('%f %f %f\n' %(v[0],v[1],v[2]))
                
        #
        # Then data related to bonds
        # 
        f.write( 'CELL_DATA %i\n' %self.N_bonds )
        f.write('SCALARS bonds double 1\nLOOKUP_TABLE default\n')       
        for bond in self.bonds:
            f.write( '%f\n' %bond[2] )
            
            
        #
        # For atom colors
        #
        #f.write('\nPOINT DATA 2 \nSCALARS colors float 1\n')
        #f.write('LOOKUP_TABLE atom_color_table\n')
        #f.write('1.0\n6.0\nLOOKUP_TABLE atom_color_table 2\n')
        #f.write('0.1 0.1 0.1 1.0\n')
        #f.write('0.7 0.7 0.7 1.0\n')
        
        f.close()  
       
if __name__=='__main__':
    import sys
    el=Element('Au')
    el.read_element_info('elements.dat')
    el()
    
    sys.exit()
    
    mol=Molecule()
    mol.input_xyz('end2.xyz')
    g= mol.pair_distr_function(0.0,20.0)
    print mol.av_bond_length()
    print mol.nr_bonds()
    f=open('test','w')
    f.writelines([str(x)+' '+str(y)+'\n' for x, y in g])
    f.close()

        
            
            



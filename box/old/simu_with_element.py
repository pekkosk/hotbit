import box.mix as mix
import numpy as nu
from box.mix import error_exit as err
from mix import type_check
from box.mix import todo
vector=nu.array

chemical_elements={}
chemical_elements['H']={'symbol':'H',\
               'name':'hydrogen',\
               'Z':1,'mass':1.008,\
               'R_cov':0.705,\
               'R_vdw':2.268,\
               'EA':0.0277,\
               'IE':0.4995}
chemical_elements['C']={'symbol':'C',\
               'name':'carbon',\
               'Z':6,'mass':12.01,\
               'R_cov':1.46,\
               'R_vdw':3.213,\
               'EA':0.0586,\
               'IE':0.4137}
chemical_elements['Au']={'symbol':'Au',\
                'name':'gold',\
                'Z':79,\
                'mass':196.970,\
                'R_cov':2.72,\
                'R_vdw':3.137,\
                'EA':0.0838,\
                'IE':0.3389}


class Element:
    """
    Class for elements. 
    
    Contains element-specific information such
    as mass, name, ionization potentials etc.
    """
    def __init__(self,symbol):
        """ Initialize element with given symbol. """
        self.data=chemical_elements[symbol.strip()]
               
    def __call__(self):
        """ Print the element symbol. """
        return self.data['symbol']
    
    def get_property(self,name):
        if name in self.data:
            return self.data[name]
        else:
            return None
        
    def pretty_print(self):
        """ Pretty output of the element properties. """
        lines=[]
        for data in self.data:
            lines.append( '%s = %s' %(data,str(self.data(data))) )
        
    def get_symbol(self): 
        return self.data['symbol']
    
    def get_mass(self):
        return self.get_property('mass')
    
    def list_integer_properties(self):
        l=[]
        for property in self.data:
            if type(self.data[property])==type(int):
                l.append(property)
        return l
    
    def list_real_properties(self):
        l=[]
        for property in self.data:
            if type(self.data[property])==type(float):
                l.append(property)
        return l
    
    
   
class Atom:
    """
    Class for atoms. 
    
    Atom consists the element information, using the class
    'Element', but contains also position, velocity etc. variable information.
    """
    def __init__(self,elem,r=vector([0.,0.,0.]),v=vector([0.,0.,0.]),\
                            f=vector([0.,0.,0.])):
        """
        Initialize atom information, element symbol, position velocity, force
        """
        if elem not in chemical_elements: 
            err('element '+elem+' not listed')
        self.element=Element(elem.strip())
        self.vectors={'r':r,'v':v,'f':f}
        self.reals={}
        self.ints={'flag':1}
               
    def __call__(self):
        """ Print some element info. """
        print 'element', self.element.get_symbol()
        for i in self.ints:      
            print i, self.ints[i]
        for r in self.reals:     
            print r, self.reals[i]
        for v in self.vectors: 
            print v, self.vectors[v]
               
    #
    # Get-methods
    #
    def get_symbol(self):
        return self.element.get_symbol()
    
    def get_element(self):
        return self.element
    
    def get_position(self):
        return self.vectors['r']
    
    def get_force(self,f):
        return self.get_vector('f')
    
    def get_velocity(self,v):
        return self.get_vector('v')
    
    def get_vector(self,name):
        return self.vectors[name]
    
    def get_int(self,name):
        if name in self.ints:
            return self.ints[name]
        else:
            return self.element.get_property(name)
    
    def get_real(self,name):
        if name in self.reals:
            return self.reals[name]
        else:
            return self.element.get_property(name)
    
    def get_mass(self):
        return self.element.get_mass()
    #
    # Set-methods
    #
    def set_vector(self,name,vec):
        type_check(name,'str')
        type_check(vec,'vector')
        self.vectors[name]=vec
    
    def set_int(self,name,i):
        type_check(i,'int')
        self.ints[name]=i
        
    def set_real(self,name,r):
        type_check(r,'float')
        assert type(r)==type(float)
        self.reals[name]=r
        
    def add_force(self,f):
        type_check(f,'vector')
        self.reals['f']=self.reals['f']+f
        
    def set_position(self,r):
        type_check(r,'vector')
        self.set_vector('r',r)
    
    def list_real_properties(self):
        return list(self.reals)+self.element.list_real_properties()
    
    def list_integer_properties(self):
        return list(self.ints)+self.element.list_integer_properties()
    
    def list_vector_properties(self):
        return list(self.vectors)
    

    
class Molecule:
    """ Class for molecules. """
    
    def __init__(self,atoms=[]):
        """ Initialize molecule. """
        self.ints={}
        self.reals={'charge':0.0,'energy':None}
        self.vectors={}
        self.tensors={}
        self.atoms=atoms
        self.N=len(atoms)
        self.N_bonds=0
        self.bonds=[]
        self.cell=[]
        self.periodic=[]
        for i in range(self.N):
            self.atoms[i].set_int('index',i)
            
    def __call__(self):
        """ Pretty print of some properties. """
        print 'N=',self.N
        for i in self.ints:      print i, self.ints[i]
        for r in self.reals:     print r, self.reals[i]
        for vec in self.vectors: print vec, self.vectors[vec]    
            
    def __add__(self, other):
        """ 
        Add (one) atom or another molecule. 
        Note: completely new molecule constructed.
        """
        new=self.atoms
        if type(other)==type(Molecule()):
            for atom in other:
                new.append(atom)
        else:
            new.append(other)
        return Molecule(new)
     
    def __iter__(self):
        """ Iterator over atoms. """
        self.itr+=1
        if self.itr==self.N: 
            raise StopIteration
        self.itr+=1
        return self.atoms[self.itr]
        
    def __iter__(self):
        """ Initialize iterator. """
        self.itr=-1
        return self
        
    #
    # Get-methods
    #
    def get_N(self): 
        return self.N
    
    def get_symbols(self):
        return [atom.get_symbol() for atom in self.atoms]
    
    def get_atom(self,i):
        return self.atoms[i]
        
    def get_vectors(self,name,mode='default'):
        v=[]
        for atom in self.atoms:
            v.append( atom.get_vector(name) )
        return mix.select_vector_mode(v,mode)

    def get_positions(self,mode='default'):
        return self.get_vectors('r',mode)
        
    def get_forces(self,mode='default'):
        return self.get_vectors('f',mode)

    def get_energy(self):
        return self.reals['energy']
   
    def get_charge(self):
        return self.reals['charge']
    #
    # Set-methods
    #
    def set_cell(vectors=None,periodic=None):
        if vectors!=None:
            self.cell=vectors
        if periodic!=None:
            self.periodic=periodic
        
    def set_real(self,key,value):
        type_check(value,'float')
        self.reals[key]=value
   
    def set_energy(self,e):   
        """ Set energy. Check also that it is 'reasonable' (approx. in any unit system). """
        if e<-100000 or 100000<e:
            err('The value for energy is suspicious.')
        self.set_real('energy',e)
        
    def set_charge(self,q):
        """ Set charge. Check also that it is reasonable. """
        if q<-100 or 100<q:
            err('The value for charge is suspicious.')
        self.set_real('charge',q)
            
    def set_vectors(self,name,v,mode='default'):
        v2=mix.select_vector_mode(v,mode)
        for atom,av in zip(self.atoms,v2):
            atom.set_vector(name,av)       
            
    def set_positions(self,r,mode='default'):
        self.set_vectors('r',r,mode)
                    
    def set_velocities(self,v,mode='default'):
        self.set_vectors('v',v,mode)
                
    def set_forces(self,f,mode='default'):
        self.set_vectors('f',f,mode)
            
    def add_forces(self,f,mode='default'):
        f2=mix.select_vector_mode(f,mode)
        for atom,fa in zip(self.atoms,f2):
            atom.add_forces(fa)            
            
    
    #
    # More complicated function
    #   
    def get_coordination_numbers(self):
        todo('def get_coordination_numbers')

    def moments_of_inertia(self):
        todo('moments_of_inertia')

    def rotate_about_axis(self,axis,angle):
        todo('rotate_about_axis')
    
    def kinetic_energy(self):
        todo('kinetic_energy')
        
    def total_mass(self):
        m_tot=0.0
        for atom in self.atoms:
            m_tot+=atom.get_mass()
        return m_tot
        
    def center_of_mass(self):
        r=vector([0.,0.,0.])
        m_tot=self.total_mass()
        for atom in self.atoms:
            r=r+atom.get_mass()*atom.get_position()
        return r/m_tot
        
    def total_momentum(self):
        todo('total_momentum')
    
    def dipole_moment(self):
        todo('dipole_moment')
    
    def translate(self):
        todo('translate')
        
    def scale_coordinates(self,x):
        todo('scale')
            
    def distance(self,i,j):
        """ 
        Return distance between points or atoms (class Atom) i
        and j (in possibly periodic cell) 
        """
        rij=self.r_ij(i,j)
        return nu.sqrt( rij[0]**2+rij[1]**2+rij[2]**2 )
    
    def r_ij(self,i,j):
        """ 
        Vector j-i with i and j positions or atoms (class Atom), 
        accounting for periodic images in minimum image convention. 
        """
        assert type(i)==type(j)
        if type(i)==type(vector([])):
            ri=i
            rj=j
        else:
            ri=i.get_position()
            rj=j.get_position()
        return rj-ri #modify this for periodic stuff ...
           
    def construct_bonds(self,factor=1.4):
        """ 
        Make the bonding list for the molecule.
        Use estimates for bond lengths from van der Waals radii.
        Make bond if R<( R_cov(1)+R_cov(2) )*factor
        """
        nb=0
        self.bonds=[]
        for i in range(self.N):
            for j in range(i+1,self.N):
                ai = self.atoms[i]
                aj = self.atoms[j]
                ri = ai.get_element().get_property('R_cov')
                rj = aj.get_element().get_property('R_cov')
                R  = self.distance( self.atoms[i],self.atoms[j] )
                if R<(ri+rj)*factor: 
                    self.bonds.append( {'i':i,'j':j,'length':R} )
                    nb+=1
        self.N_bonds = nb
           
    #def set_bonds(self,file):
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

    def write_vtk(self,file):
        """ vtk output of molecule (all scalar and vector properties) """
        f=open(file,'w')
        f.write('# vtk DataFile Version 2.0 \nNOTB simulation\n')
        f.write('ASCII \nDATASET UNSTRUCTURED_GRID\n')
        f.write('POINTS %i double \n ' %self.N)
        
        # Point data (atom coordinates) and cell data (bonds)
        if self.N_bonds==0: 
            self.construct_bonds()
        for atom in self.atoms:
            s=a2s( atom.get_position() )
            f.write('%s\n' %s )
        f.write( 'CELLS %i %i\n' %(self.N_bonds,3*self.N_bonds) )
        for bond in self.bonds:
            f.write( '2 %i %i\n' %(bond['i'],bond['j']) )
        f.write( 'CELL_TYPES %i\n' %self.N_bonds )
        for bond in self.bonds:
            f.write( '3\n' )    
           
        # First the data related to atoms
        f.write('POINT_DATA %i\n' %self.N)
        
        real_properties=self.atoms[0].list_real_properties()
        for real in real_properties:
            f.write('SCALARS %s double 1\nLOOKUP_TABLE default\n' %real)       
            for atom in self.atoms:
                f.write('%f\n' %atom.get_real(real))
        
        int_properties=self.atoms[0].list_integer_properties()
        for i in int_properties:
            f.write('SCALARS %s double 1\nLOOKUP_TABLE default\n' %i)       
            for atom in self.atoms:
                f.write('%f\n' %atom.get_int(i))
        
        vec_properties=self.atoms[0].list_vector_properties()
        for vec in vec_properties:
            f.write('SCALARS %s double 1\nLOOKUP_TABLE default\n' %vec)       
            for atom in self.atoms:
                v=atom.get_vector(vec)
                f.write( '%f %f %f\n' %(v[0],v[1],v[2]) )
        
        # Then data related to bonds
        f.write( 'CELL_DATA %i\n' %self.N_bonds )
        f.write('SCALARS bonds double 1\nLOOKUP_TABLE default\n')       
        for bond in self.bonds:
            f.write( '%f\n' %bond['length'] )
            
        f.close()  
       
    
if __name__=='__main__':
    x1=vector([1.,0.,0.])
    x2=vector([2.,0.,0.])
    x3=vector([3.,0.,0.])
    a1=Atom('H',x1)
    a2=Atom('C',x2)
    a3=Atom('C',x3)
    mol=Molecule([a1,a2,a3])
    print mol.get_symbols()
    r_list=[x1,x2,x3]
    mol.set_positions(r_list)
    print mol.get_positions()
    print mol.total_mass()
    print mol.center_of_mass()
    print mol.distance(x1,x2)
    print mol.distance(a1,a2)
    #print mol.pair_distribution_list()
    #print mol.pair_distribution_function(-1.0,4.0)
    #print 'average bond length',mol.average_bond_length()
    #mol.write_vtk('koe.vtk')
    #print mol.center_of_mass()

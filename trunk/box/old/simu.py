import box.mix as mix
import numpy as nu
from box.mix import error_exit as err
from mix import type_check
from box.mix import todo
from copy import copy
vector=nu.array

chemical_elements={}
chemical_elements['H']={'symbol':'H',\
               'name':'hydrogen',\
               'Z':1,\
               'mass':1.008,\
               'R_cov':0.705,\
               'R_vdw':2.268,\
               'EA':0.0277,\
               'IE':0.4995}
chemical_elements['C']={'symbol':'C',\
               'name':'carbon',\
               'Z':6,\
               'mass':12.01,\
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


#class Element:
    #"""
    #Class for elements. 
    
    #Contains element-specific information such
    #as mass, name, ionization potentials etc.
    #"""
    #def __init__(self,symbol):
        #""" Initialize element with given symbol. """
        #self.data=chemical_elements[symbol.strip()]
               
    #def __call__(self):
        #""" Print the element symbol. """
        #return self.data['symbol']
    
    #def get_property(self,name):
        #if name in self.data:
            #return self.data[name]
        #else:
            #return None
        
    #def pretty_print(self):
        #""" Pretty output of the element properties. """
        #lines=[]
        #for data in self.data:
            #lines.append( '%s = %s' %(data,str(self.data(data))) )
        
    #def get_symbol(self): 
        #return self.data['symbol']
    
    #def get_mass(self):
        #return self.get_property('mass')
    
    #def list_integer_properties(self):
        #l=[]
        #for property in self.data:
            #if type(self.data[property])==type(int):
                #l.append(property)
        #return l
    
    #def list_real_properties(self):
        #l=[]
        #for property in self.data:
            #if type(self.data[property])==type(float):
                #l.append(property)
        #return l
    
    
   
class Atom:
    """
    Class for atoms. 
    
    Atom consists the element information, using the class
    'Element', but contains also position, velocity etc. variable information.
    """
    def __init__(self,symbol,r=vector([0.,0.,0.]),v=vector([0.,0.,0.]),\
                            f=vector([0.,0.,0.])):
        """
        Initialize atom information, element symbol, position velocity, force
        """
        if symbol not in chemical_elements: 
            err('element '+symbol+' not listed')
        
        self.data=copy( chemical_elements[symbol.strip()] )
        self.data['r']=r
        self.data['v']=v
        self.data['f']=f
        self.data['flag']=1
               
    def __call__(self):
        """ Print some element info. """
        print 'element', self.data['symbol']
        for d in self.data:
            print d,'=',seld.data[d]
        
    #
    # Get-methods
    #
    def get_property(self,name):
        return self.data[name]
    
    def get_symbol(self):
        return self.get_property('symbol')
    
    def get_position(self):
        return self.get_property('r')
    
    def get_force(self):
        return self.get_property('f')
    
    def get_velocity(self):
        return self.get_property('v')
        
    def get_mass(self):
        return self.get_property('mass')
    #
    # Set-methods
    #
    def set_property(self,name,p):
        self.data[name]=p
                
    def set_position(self,r):
        type_check(r,'vector')
        self.set_property('r',r)
        
    def set_velocity(self,v):
        type_check(v,'vector')
        if mix.norm(v)>1000:
            self()
            err('Force suspiciously large.')
        self.set_property('v',v)
        
    def set_force(self,f):
        type_check(f,'vector')
        if mix.norm(f)>1000: 
            self()
            err('Force suspiciously large.')
        self.set_property('f',f)
        
    def add_force(self,f):
        type_check(f,'vector')
        set_force( self.get_force()+f )
        
    def list_properties(self):
        return list(self.data)
    

    
class Molecule:
    """ Class for molecules. """
    
    def __init__(self,atoms=[]):
        """ Initialize molecule. """
        self.data={'charge':0.0,'energy':None}
        self.atoms=atoms
        self.N=len(atoms)
        self.N_bonds=0
        self.bonds=[]
        self.cell=[]
        self.periodic=[]
        for i in range(self.N):
            self.atoms[i].set_property('index',i)
            
    def __call__(self):
        """ Pretty print of some properties. """
        print 'N=',self.N
        for d in self.data:
            print d,'=', self.data[d]
            
    def __add__(self, other):
        """ 
        Add (one) atom or another molecule. 
        Note: completely new molecule constructed.
        """
        new=copy(self.atoms)
        if type(other)==type(Molecule()):
            for atom in other:
                new.append(atom)
        else:
            new.append(other)
        return Molecule(new)
     
    def __iter__(self):
        """ Initialize iterator. """
        self.itr=-1
        return self 
     
    def next(self):
        """ Iterator over atoms. """
        self.itr+=1
        if self.itr==self.N: 
            raise StopIteration
        return self.atoms[self.itr]
        
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
            va=atom.get_property(name)
            assert type(va)==type(vector([]))
            v.append(va)
        return mix.select_vector_mode(v,mode)

    def get_positions(self,mode='default'):
        return self.get_vectors('r',mode)
        
    def get_forces(self,mode='default'):
        return self.get_vectors('f',mode)

    def get_velocities(self,mode='default'):
        return self.get_vectors('v',mode)
    
    def get_energy(self):
        return self.data['energy']
   
    def get_charge(self):
        return self.data['charge']
    #
    # Set-methods
    #
    def set_cell(vectors=None,periodic=[False,False,False]):
        """ Set periodic cell with given (three) vectors and periodicities in 3 directions. """
        if vectors!=None:
            self.cell=vectors
        if periodic!=None:
            self.periodic=periodic
        
    def set_property(self,name,value):
        self.data[name]=value
       
    def set_energy(self,e):   
        """ Set energy. Check also that it is 'reasonable' (approx. in any unit system). """
        if e<-100000 or 100000<e:
            err('The value for energy is suspicious.')
        self.set_property('energy',e)
        
    def set_charge(self,q):
        """ Set charge. Check also that it is reasonable. """
        if q<-100 or 100<q:
            err('The value for charge is suspicious.')
        self.set_property('charge',q)
            
    def set_vectors(self,name,v):
        v2=mix.select_vector_mode(v,mode='default')
        for atom,av in zip(self.atoms,v2):
            atom.set_property(name,av)       
            
    def set_positions(self,r):
        self.set_vectors('r',r)
                    
    def set_velocities(self,v):
        self.set_vectors('v',v)
                
    def set_forces(self,f):
        self.set_vectors('f',f)
            
    def add_forces(self,f):
        f2=mix.select_vector_mode(f,mode='default')
        for atom,fa in zip(self.atoms,f2):
            atom.add_force(fa)            
            
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
                ri = ai.get_property('R_cov')
                rj = aj.get_property('R_cov')
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
            s=mix.a2s( atom.get_position() )
            f.write('%s\n' %s )
        f.write( 'CELLS %i %i\n' %(self.N_bonds,3*self.N_bonds) )
        for bond in self.bonds:
            f.write( '2 %i %i\n' %(bond['i'],bond['j']) )
        f.write( 'CELL_TYPES %i\n' %self.N_bonds )
        for bond in self.bonds:
            f.write( '3\n' )    
           
        # First the data related to atoms
        f.write('POINT_DATA %i\n' %self.N)
        for property in self.atoms[0].list_properties():
            tp=type( self.atoms[0].get_property(property) )
            if tp==type(1) or tp==type(1.0):
                print>>f, 'SCALARS %s double 1\nLOOKUP_TABLE default' %property
                for atom in self.atoms:
                    print>>f, atom.get_property(property)
            elif tp==type(vector([])):
                print>>f, 'VECTORS %s double\nLOOKUP_TABLE default' %property
                for atom in self.atoms:
                    print>>f, mix.a2s(atom.get_property(property))
        
        # Then data related to bonds
        print>>f, 'CELL_DATA %i' %self.N_bonds 
        print>>f, 'SCALARS bond_length double 1\nLOOKUP_TABLE default'
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

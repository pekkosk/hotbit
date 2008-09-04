import os
import numpy as nu
from box import mix
from box.interpolation import Function
import numpy as npy
from box.data import data
from box.data import atom_valence
from copy import copy
find=mix.find_value
                
class Element:
    def __init__(self,element):
        """
        Initialize element with given symbol or from file.
        ('H' or 'H.elm'; .elm file first searched from present directory,
        then from the default HOTBIT_PARAMETERS/element.elm
        """
        i=element.find('elm')

        if i>=0:
            # element is full file name; get symbol
            filename=element
            j=element.rfind('/')
            if j>=0:
                ind=j+1
            else:
                ind=0
            symbol=element[ind]
            if element[ind+1].islower() and element[ind+1].isalpha():
                symbol+=element[ind+1]
        else:
            symbol=element
            filename=element+'.elm'
            
        self.data=copy(data[symbol])
        
        # if .elm file exists, read more data
        default='%s/%s' %(os.environ.get('HOTBIT_PARAMETERS'),filename)
        if os.path.isfile(filename):
            file=filename
        else:
            file=default
            
        self.data['file']=file
        f=open(file,'r')
        symb=find(f,'symbol').strip()
        assert symb==symbol
        self.data['valence']=float(find(f,'valence'))
        self.data['orbitals']=int(find(f,'orbitals'))
        self.data['U']=float(find(f,'U'))
        self.data['FWHM']=float(find(f,'FWHM'))
        energies=find(f,'energies',fmt='onestr').split()
        self.data['energies']=npy.array([float(e) for e in energies])
        self.data['comment']=find(f,'comment',fmt='onestr',default='')
        self.lst=['s','px','py','pz','dxy','dyz','dzx','dx2-y2','d3z2-r2']
        self.valence_orbitals=copy(atom_valence[symbol])
        self.read_functions(file)
        
        
    def set(self,key,value):
        self.data[key]=value
                
    def get_comment(self):
        """ Get comments describing the present element. """
        return self.data['comment']                        
    
    def _transform_orbital(self,orb):
        """ Transform 'px->'2p' etc. """
        if orb in self.valence_orbitals:
            return orb
        for v in self.valence_orbitals:
            if orb[0]==v[1]: return v    
        
    def get_valence_orbitals(self):
        """ Return ['2s','2p']... """
        return self.valence_orbitals
            
    def get_file(self):
        return self.data['file']
            
    def get_symbol(self):
        return self.data['symbol']
    
    def get_valence(self):
        return self.data['valence']
   
    def get_energies(self):
        return self.data['energies']
    
    def get_epsilon(self,nl):
        return self.data['epsilon'][nl]
        
    def get_comments(self):
        return self.data['comments']
    
    def get_nr_orbitals(self):
        return self.data['orbitals']
    
    def get_orbitals(self):
        no=self.get_nr_orbitals()
        return self.lst[:no]
    
    def get_U(self):
        return self.data['U']
   
    def get_FWHM(self):
        return self.data['FWHM']
   
    def get_Z(self):
        return self.data['Z']
        
    def write_to_file(self,file):
        """
        Write simple element file (.elm)
        """
        f=open(file,'w')
        for item in self.data:
            print>>f, item,'=',self.data[item]
        f.close()
        
    def unl(self,r,nl,der=0):
        """ Return unl=Rnl*r. """
        nl=self._transform_orbital(nl)
        return self.functions['unl'][nl](r,der=der)
    
    def Rnl(self,r,nl,der=0):
        """ Return Rnl(r). """
        nl=self._transform_orbital(nl)
        return self.functions['Rnl'][nl](r,der=der)
    
    def get_Rnl(self,nl):
        """ Return the actual Rnl function object. """
        nl=self._transform_orbital(nl)
        return self.functions['Rnl'][nl]  
    
    def v_effective(self,r):
        return self.functions['v_effective'](r)
    
    def confinement(self,r):
        return self.functions['confinement'](r)
        
    def read_functions(self,file):
        """ 
        Read radial wave functions (R=u/r), self-consistent potentials, 
        confinements, etc. from given file. 
        """
        o=open(file)        
        self.functions={'unl':{},'Rnl':{}}
        for nl in self.get_valence_orbitals():
            m=find(o,'u_%s' %nl,fmt='matrix',default=nu.array([[0,0],[1,0],[2,0],[3,0]]))  
            self.functions['unl'][nl]=Function('spline',m[:,0],m[:,1])
            self.functions['Rnl'][nl]=Function('spline',m[:,0],m[:,1]/m[:,0])                                    
            
        #m=find(o,'v_effective',fmt='matrix')  
        #self.functions['v_effective']=Function('spline',m[:,0],m[:,1])
        #m=find(o,'confinement',fmt='matrix')  
        #self.functions['confinement']=Function('spline',m[:,0],m[:,1])
        
        
        
        
if __name__=='__main__':
    elm=Element('H.elm')
    print elm.get_symbol()
    print elm.__dict__
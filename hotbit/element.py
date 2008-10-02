import os
import numpy as nu
from box import mix
from box.interpolation import Function
import numpy as npy
from box.data import data
from copy import copy
find=mix.find_value

orbital_list=['s','px','py','pz','dxy','dyz','dzx','dx2-y2','d3z2-r2']
                
class Element:
    def __init__(self,element):
        """
        Initialize element with given symbol (e.g. 'C') or from file (e.g. 'C.elm').
        (elm file first searched from present directory,
        then from the default HOTBIT_PARAMETERS/[symbol].elm
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
        
        # read data from elm-file
        self.data['valence_number']=copy( data[symbol]['valence_number'] )
        self.data['valence_orbitals']=copy( data[symbol]['valence_orbitals'] )
        self.data['U']=float(find(f,'U'))
        self.data['FWHM']=float(find(f,'FWHM'))
        self.data['epsilon']={}
        for orbital in self.data['valence_orbitals']:
            self.data['epsilon'][orbital]=float( find(f,'epsilon_%s' %orbital) )
        self.data['comment']=find(f,'comment',fmt='lines',default='no comment')
        self._read_functions(file)
        
        # for internal...
        energies=[]            
        for nl in self.data['valence_orbitals']:            
            eps=self.get_epsilon(nl)
            if   's' in nl: n=1
            elif 'p' in nl: n=3
            elif 'd' in nl: n=5
            energies.extend( [eps]*n )                
        self.data['onsite_energies']=energies                
        self.data['nr_basis_orbitals']=len(energies)                
        self.data['valence_energies']=npy.array([float(e) for e in energies])
        
        
    def set(self,key,value):
        self.data[key]=value
                
                
    def get_comment(self):
        """ Get comments describing the present element. """
        return self.data['comment']                        
    
    
    def _transform_orbital(self,orb):
        """ Transform 'px->'2p' etc. """
        if orb in self.data['valence_orbitals']:
            return orb
        for v in self.data['valence_orbitals']:
            if orb[0]==v[1]: return v    
        
        
    def get_valence_orbitals(self):
        """ Return ['2s','2p']... """
        return self.data['valence_orbitals']
            
            
    def get_symbol(self):
        """ Return 'H', 'C', etc. """
        return self.data['symbol']
    
    
    def get_valence_number(self):
        """ Return the number of valence electrons. """
        return self.data['valence_number']
   
   
    def get_epsilon(self,nl):
        """ Return on-site energy for given orbital '2p',... """
        return self.data['epsilon'][nl]
        
        
    def get_comments(self):
        """ Return comments regarding the element file. """
        return self.data['comment']
    
    
    def get_nr_basis_orbitals(self):
        """ Return the total number of basis orbitals (1,4 or 9) """
        return self.data['nr_basis_orbitals']
    
    
    def get_onsite_energies(self):
        """ Return on-site energies for all basis functions. """
        return self.data['onsite_energies']
    
    
    def get_orbital_types(self):
        """ Return list of valence orbital types ['s','px','py',...] """
        no=self.get_nr_basis_orbitals()
        return orbital_list[:no]
    
    
    def get_U(self):
        return self.data['U']
   
   
    def get_FWHM(self):
        return self.data['FWHM']
   
   
    def get_Z(self):
        return self.data['Z']
      
      
    def get_file(self):
        """ Return file name where data was read. """
        return self.data['file']
              
        
    def unl(self,r,nl,der=0):
        """ Return unl=Rnl*r. """
        nl=self._transform_orbital(nl)
        return self.functions['unl'][nl](r,der=der)
    
    
    def Rnl(self,r,nl,der=0):
        """ Return Rnl(r). """
        nl=self._transform_orbital(nl)
        return self.functions['Rnl'][nl](r,der=der)
    
    
    def get_Rnl_function(self,nl):
        """ Return the actual Rnl function object. """
        nl=self._transform_orbital(nl)
        return self.functions['Rnl'][nl]  
    
    
    def v_effective(self,r):
        """ Return Kohn-Sham effective potential. """
        return self.functions['v_effective'](r)
    
    
    def confinement(self,r):
        """ Return the confinement potential used to generate pseudo-atom. """
        return self.functions['confinement'](r)
        
        
    def _read_functions(self,file):
        """ 
        Read radial wave functions (R=u/r), self-consistent potentials, 
        confinements, etc. from given file. 
        """
        o=open(file)        
        default=nu.array([[0,0],[1,0],[2,0],[3,0]])
        self.functions={'unl':{},'Rnl':{}}
        for nl in self.get_valence_orbitals():
            m=find(o,'u_%s' %nl,fmt='matrix',default=default)  
            self.functions['unl'][nl]=Function('spline',m[:,0],m[:,1])
            self.functions['Rnl'][nl]=Function('spline',m[:,0],m[:,1]/m[:,0])                                    
            
        m=find(o,'v_effective',fmt='matrix',default=default)  
        self.functions['v_effective']=Function('spline',m[:,0],m[:,1])
        m=find(o,'confinement',fmt='matrix',default=default)  
        self.functions['confinement']=Function('spline',m[:,0],m[:,1])
        o.close()
        
    def write_to_file(self,file):
        """
        Write simple element file (.elm) 
        """
        f=open(file,'w')
        for item in self.data:
            print>>f, item,'=',self.data[item]
        f.close()        
        
        
if __name__=='__main__':
    elm=Element('H.elm')
    print elm.get_symbol()
    print elm.__dict__
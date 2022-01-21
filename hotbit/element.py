from __future__ import print_function

import os
import numpy as np
from box import mix
import numpy as npy
from box.data import data
from hotbit.io import read_element

orbital_list=['s','px','py','pz','dxy','dyz','dzx','dx2-y2','d3z2-r2']
                
class Element:
    def __init__(self,element):
        """
        Initialize element with given symbol (e.g. 'C') or from file (e.g. 'C.elm').
        (elm file first searched from present directory,
        then from the default HOTBIT_PARAMETERS/[symbol].elm
        """
        i=element.rfind('.')+1
        
        if i>0:
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
            
        # if .elm file exists, read more data
        hbp = os.environ.get('HOTBIT_PARAMETERS')
        default='%s/%s' %(hbp,filename)
        if os.path.isfile(filename):
            file=filename
        elif os.path.isfile(default):
            file=default
        else:
            raise AssertionError('Element data for %s not found neither in %s '
                                 'nor in %s' % (element,'.',hbp))

        self.data, self.functions = read_element(file, symbol)
        self.data['file']=file

        
    def set(self,key,value):
        self.data[key]=value
        
        
    def get_free_atom_energy(self):
        '''
        Return the energy of a free element (sum_i f_i*e_i)
        '''
        e=0.0
        for valence in self.data['valence_orbitals']:
            e += self.data['configuration'][valence]*self.get_epsilon(valence)
        return e
                
                
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


    def get_C6(self):
        return float(self.data['C6'])


    def get_p(self):
        return float(self.data['p'])


    def get_R0(self):
        return float(self.data['R0'])


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


    def has_Rnl_functions(self):
        """ Does this element has basis function information? """
        return self.functions is not None
    
    
    def effective_potential(self,r):
        """ Return Kohn-Sham effective potential. """
        return self.functions['effective_potential'](r)
    
    
    def confinement_potential(self,r):
        """ Return the confinement potential used to generate pseudo-atom. """
        return self.functions['confinement_potential'](r)
        
        
    def get_wf_range(self,nl,fractional_limit=1E-7):
        """ Return the maximum r for which |R(r)|<fractional_limit*max(|R(r)|) """
        rgrid=np.linspace(1,20,200)
        rnl=np.array([self.Rnl(r,nl) for r in rgrid])
        wfmax=max(abs(rnl))
        for r,wf in zip(rgrid[-1::-1],rnl[-1::-1]):
            if abs(wf)>fractional_limit*wfmax: 
                return r
                
        
    def write_to_file(self,file):
        """
        Write simple element file (.elm) 
        """
        f=open(file,'w')
        for item in self.data:
            print(item,'=',self.data[item], file=f)
        f.close()


    def update_vdw(self,p,R0):
        self.data['p'] = p
        self.data['R0'] = R0
        self.data['C6'] = None

        
if __name__=='__main__':
    elm=Element('H.elm')
    print(elm.get_symbol())
    print(elm.__dict__)

# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

data={}

data['H'] ={'symbol':'H',\
            'name':'hydrogen',\
            'Z':1,\
            'mass':1.00797,\
            'R_cov':0.37,\
            'R_vdw':1.2}
            
data['C']={'symbol':'C',\
            'name':'carbon',\
            'Z':6,\
            'mass':12.01115,\
            'R_cov':0.77,\
            'R_vdw':1.7}
            
data['Au']={'symbol':'Au',\
            'name':'gold',\
            'Z':79,\
            'mass':196.967,\
            'R_cov':1.44,\
            'R_vdw':1.66}
            
data['N'] ={'symbol':'N',\
            'name':'nitrogen',\
            'Z':7,\
            'mass':14.0067,\
            'R_cov':0.75,\
            'R_vdw':1.55}
            
data['O'] ={'symbol':'O',\
            'name':'oxygen',\
            'Z':8,\
            'mass':15.9994,\
            'R_cov':0.73,\
            'R_vdw':1.52}
            
data['Na']={'symbol':'Na',\
            'name':'sodium',\
            'Z':11,\
            'mass':22.9898,\
            'R_cov':1.54,\
            'R_vdw':2.27}
           
data['K']={'symbol':'K',\
            'name':'potassium',\
            'Z':19,\
            'mass':39.102,\
            'R_cov':2.03,\
            'R_vdw':2.75} 
            
data['Ti']={'symbol':'Ti',\
            'name':'titanium',\
            'Z':22,\
            'mass':47.90,\
            'R_cov':1.36,\
            'R_vdw':2.00}
            
data['Cl']={'symbol':'Cl',\
            'name':'chlorine',\
            'Z':17,\
            'mass':35.453,\
            'R_cov':1.02,\
            'R_vdw':1.75}
            
data['Mg']={'symbol':'Mg',\
            'name':'magnesium',\
            'Z':12,\
            'mass':24.305,\
            'R_cov':1.41,\
            'R_vdw':1.73}            
                        
            
atom_valence={}
atom_valence['H'] =['1s']
atom_valence['He']=['1s']
#second row
atom_valence['Li']=['2s','2p']
atom_valence['Be']=['2s','2p']
atom_valence['B'] =['2s','2p']
atom_valence['C'] =['2s','2p']
atom_valence['N'] =['2s','2p']
atom_valence['O'] =['2s','2p']
atom_valence['F'] =['2s','2p']
atom_valence['Ne']=['2s','2p']
#third row
atom_valence['Na']=['3s','3p']
atom_valence['Mg']=['3s','3p']
atom_valence['Cl']=['3s','3p']
# fourth row
atom_valence['Ti']=['3d','4s','4p']
atom_valence['Au']=['6s','6p','5d']


aux=[ ['H', '',{'1s':1}],\
      ['He','',{'1s':2}],\
      # second row
      ['Li','He',{'2s':1,'2p':0}],\
      ['Be','He',{'2s':2,'2p':0}],\
      ['B', 'He',{'2s':2,'2p':1}],\
      ['C', 'He',{'2s':2,'2p':2}],\
      ['N', 'He',{'2s':2,'2p':3}],\
      ['O', 'He',{'2s':2,'2p':4}],\
      ['F', 'He',{'2s':2,'2p':5}],\
      ['Ne','He',{'2s':2,'2p':6}],\
      # third row
      ['Na','Ne',{'3s':1,'3p':0}],\
      ['Mg','Ne',{'3s':2,'3p':0}],\
      ['Cl','Ne',{'3s':2,'3p':5}],\
      ['Ar','Ne',{'3s':2,'3p':6}],\
      # fourth row
      ['Ti','Ar',{'3d':2,'4s':2,'4p':0}],\
      ['Kr','Ar',{'3d':10,'4s':2,'4p':6}],\
      # fifth row
      ['Xe','Kr',{'4d':10,'5s':2,'5p':6}],\
      # sixth row
      ['Au','Xe',{'4f':14,'5d':10,'6s':1,'6p':0}] ]
          
atom_occupations={}          
for item in aux:
    el, core, occu=item
    if core!='': 
        atom_occupations[el]=atom_occupations[core].copy()
    else:
        atom_occupations[el]={}        
    atom_occupations[el].update(occu)
    

#atom_valence['']=[]
            
            #chemical_elements['H']={'symbol':'H',\
               #'name':'hydrogen',\
               #'Z':1,\
               #'mass':1.008,\
               #'R_cov':0.705,\
               #'R_vdw':2.268,\
               #'EA':0.0277,\
               #'IE':0.4995}
#chemical_elements['C']={'symbol':'C',\
               #'name':'carbon',\
               #'Z':6,\
               #'mass':12.01,\
               #'R_cov':1.46,\
               #'R_vdw':3.213,\
               #'EA':0.0586,\
               #'IE':0.4137}
#chemical_elements['Au']={'symbol':'Au',\
                #'name':'gold',\
                #'Z':79,\
                #'mass':196.970,\
                #'R_cov':2.72,\
                #'R_vdw':3.137,\
                #'EA':0.0838,\
                #'IE':0.3389}
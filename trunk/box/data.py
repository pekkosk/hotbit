# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.
#
# Experimental data * mass, R_cov (2008 data), R_vdw, EA from www.webelements.com
#                   * IE from gElemental 1.2.0
#
from numpy import nan

data={}

data['H'] ={'Z':1, 'symbol':'H', 'name':'hydrogen', 'mass': 1.0079,  'R_cov':0.31,'R_vdw':1.20,'IE':0.0135,'EA':72.27 }            
data['He']={'Z':2, 'symbol':'He'}
data['Li']={'Z':3, 'symbol':'Li'}
data['Be']={'Z':4, 'symbol':'Be'}
data['B'] ={'Z':5, 'symbol':'B'}
data['C'] ={'Z':6, 'symbol':'C', 'name':'carbon',   'mass':12.0107,  'R_cov':0.76,'R_vdw':1.70,'IE':11.256,'EA':1.594}            
data['N'] ={'Z':7, 'symbol':'N', 'name':'nitrogen', 'mass':14.0067,  'R_cov':0.71,'R_vdw':1.55,'IE':14.527,'EA':0.072}            
data['O'] ={'Z':8, 'symbol':'O', 'name':'oxygen',   'mass':15.9994,  'R_cov':0.66,'R_vdw':1.52,'IE':13.612,'EA':1.460}            
data['F'] ={'Z':9, 'symbol':'F'}
data['Ne']={'Z':10,'symbol':'Ne'}
data['Na']={'Z':11,'symbol':'Na','name':'sodium',   'mass':22.9898,  'R_cov':1.66,'R_vdw':2.27,'IE':5.136, 'EA':0.547}           
data['Mg']={'Z':12,'symbol':'Mg','name':'magnesium','mass':24.3050,  'R_cov':1.41,'R_vdw':1.73,'IE':7.642, 'EA':0.000}                        
data['S'] ={'Z':16, 'symbol':'S', 'name':'sulfur', 'mass':32.065, 'R_cov':1.02, 'R_vdw':1.80, 'IE':999.6, 'EA':200} #FIXME these are kJ/mol, what are the others?
data['Cl']={'Z':17,'symbol':'Cl','name':'chlorine', 'mass':35.4530,  'R_cov':1.02,'R_vdw':1.75,'IE':12.962,'EA':3.615}            
data['Ar']={'Z':18,'symbol':'Ar'}
data['K'] ={'Z':19,'symbol':'K', 'name':'potassium','mass':39.0983,  'R_cov':2.03,'R_vdw':2.75,'IE':4.338, 'EA':0.501}             
data['Ti']={'Z':22,'symbol':'Ti','name':'titanium', 'mass':47.8760,  'R_cov':1.60,'R_vdw':2.15,'IE':6.825, 'EA':0.078}            
data['Kr']={'Z':36,'symbol':'Kr'}
data['Xe']={'Z':54,'symbol':'Xe'}
data['Au']={'Z':79,'symbol':'Au','name':'gold',     'mass':196.9666, 'R_cov':1.36,'R_vdw':1.66,'IE':9.221, 'EA':2.308}                    
                        
# update with valence orbital data                        
valence_orbitals={}
valence_orbitals['H'] =['1s']
valence_orbitals['He']=['1s']
valence_orbitals['Li']=['2s','2p']
valence_orbitals['Be']=['2s','2p']
valence_orbitals['B'] =['2s','2p']
valence_orbitals['C'] =['2s','2p']
valence_orbitals['N'] =['2s','2p']
valence_orbitals['O'] =['2s','2p']
valence_orbitals['F'] =['2s','2p']
valence_orbitals['Ne']=['2s','2p']
valence_orbitals['Na']=['3s','3p']
valence_orbitals['Mg']=['3s','3p']
valence_orbitals['S'] =['3s','3p']
valence_orbitals['Cl']=['3s','3p']
valence_orbitals['Ar']=[]
valence_orbitals['Ti']=['3d','4s','4p']
valence_orbitals['Kr']=[]
valence_orbitals['Xe']=[]
valence_orbitals['Au']=['6s','6p','5d']
for key in valence_orbitals:
    data[key]['valence_orbitals']=valence_orbitals[key]


# Set electronic configurations (orbital occupations)
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
      ['S', 'Ne',{'3s':2,'3p':4}],\
      ['Cl','Ne',{'3s':2,'3p':5}],\
      ['Ar','Ne',{'3s':2,'3p':6}],\
      # fourth row
      ['Ti','Ar',{'3d':2,'4s':2,'4p':0}],\
      ['Kr','Ar',{'3d':10,'4s':2,'4p':6}],\
      # fifth row
      ['Xe','Kr',{'4d':10,'5s':2,'5p':6}],\
      # sixth row
      ['Au','Xe',{'4f':14,'5d':10,'6s':1,'6p':0}] ]
          
configurations={}          
for item in aux:
    el, core, occu=item
    if core!='': 
        configurations[el]=configurations[core].copy()
    else:
        configurations[el]={}        
    configurations[el].update(occu)
for key in configurations:
    config=configurations[key]    
    data[key]['configuration']=config
    data[key]['valence_number']=sum( [config[orbital] for orbital in data[key]['valence_orbitals']] )
    
if __name__=='__main__':    
    for symbol in data:
        print 'X'*40
        for d in data[symbol]:
            print d,data[symbol][d]
    
    
   
   

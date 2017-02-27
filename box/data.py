# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.
#
# Experimental data * mass, R_cov (2008 data), R_vdw, EA from www.webelements.com
#                   * IE from gElemental 1.2.0
#
# UNITS:
#     * mass in amu
#     * all radii in Angstrom
#     * all energies in eV
#
from numpy import nan

data={}

data['H'] ={'Z':1, 'symbol':'H', 'name':'hydrogen', 'mass': 1.0079,  'R_cov':0.31,'R_vdw':1.20,'IE':0.0135,'EA':72.27 }
data['He']={'Z':2, 'symbol':'He'}
data['Li']={'Z':3, 'symbol':'Li', 'R_cov':1.28}
data['Be']={'Z':4, 'symbol':'Be'}
data['B'] ={'Z':5, 'symbol':'B', 'name':'boron',    'mass':10.81,    'R_cov':0.84,'R_vdw':1E99,'IE':8.294, 'EA':0.277}
data['C'] ={'Z':6, 'symbol':'C', 'name':'carbon',   'mass':12.0107,  'R_cov':0.76,'R_vdw':1.70,'IE':11.256,'EA':1.594}
data['N'] ={'Z':7, 'symbol':'N', 'name':'nitrogen', 'mass':14.0067,  'R_cov':0.71,'R_vdw':1.55,'IE':14.527,'EA':0.072}
data['O'] ={'Z':8, 'symbol':'O', 'name':'oxygen',   'mass':15.9994,  'R_cov':0.66,'R_vdw':1.52,'IE':13.612,'EA':1.460}
data['F'] ={'Z':9, 'symbol':'F'}
data['Ne']={'Z':10,'symbol':'Ne'}
data['Na']={'Z':11,'symbol':'Na','name':'sodium',   'mass':22.9898,  'R_cov':1.66,'R_vdw':2.27,'IE':5.136, 'EA':0.547}
data['Mg']={'Z':12,'symbol':'Mg','name':'magnesium','mass':24.3050,  'R_cov':1.41,'R_vdw':1.73,'IE':7.642, 'EA':0.000}
data['Si']={'Z':14,'symbol':'Si','name':'silicon',  'mass':28.0855}
data['P'] ={'Z':15,'symbol':'P', 'name':'phosphorus','mass':30.9737, 'R_cov':1.07,'R_vdw':1.80,'IE':10.487, 'EA':0.746}
data['S'] ={'Z':16, 'symbol':'S', 'name':'sulfur',  'mass':32.065,   'R_cov':1.05,'R_vdw':1.80,'IE':10.356, 'EA':2.072}
data['Cl']={'Z':17,'symbol':'Cl','name':'chlorine', 'mass':35.4530,  'R_cov':1.02,'R_vdw':1.75,'IE':12.962,'EA':3.615}
data['Ar']={'Z':18,'symbol':'Ar'}
data['K'] ={'Z':19,'symbol':'K', 'name':'potassium','mass':39.0983,  'R_cov':2.03,'R_vdw':2.75,'IE':4.338, 'EA':0.501}
data['Ti']={'Z':22,'symbol':'Ti','name':'titanium', 'mass':47.8760,  'R_cov':1.60,'R_vdw':2.15,'IE':6.825, 'EA':0.078}
data['Kr']={'Z':36,'symbol':'Kr'}
data['Sr']={'Z':38,'symbol':'Sr','name':'strontium','mass':87.62,    'R_cov':1.95,'R_vdw':2.49, 'IE':5.69,'EA':0.052}
data['Mo']={'Z':42,'symbol':'Mo','name':'molybdenum','mass':95.94,   'R_cov':1.57,'R_vdw':2.10, 'IE':7.08,'EA':0.744}
data['Pd']={'Z':46,'symbol':'Pd'}
data['Sn']={'Z':50, 'symbol':'Sn', 'R_cov':1.39}
data['Xe']={'Z':54,'symbol':'Xe'}
data['Pt']={'Z':78,'symbol':'Pt','name':'platinum', 'mass':195.084,  'R_cov':1.36,'R_vdw':1.75,'IE':9.013, 'EA':2.127}
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
valence_orbitals['Si']=['3s','3p']
valence_orbitals['P'] =['3s','3p']
valence_orbitals['S'] =['3s','3p']
valence_orbitals['Cl']=['3s','3p']
valence_orbitals['Ar']=[]
valence_orbitals['K']=['4s','4p']
valence_orbitals['Ti']=['3d','4s','4p']
valence_orbitals['Kr']=[]
valence_orbitals['Sr']=['5s','5p','4d']
valence_orbitals['Mo']=['5s','5p','4d']
valence_orbitals['Pd']=['5s','5p','4d']
valence_orbitals['Sn']=['5s','5p']
valence_orbitals['Xe']=[]
valence_orbitals['Pt']=['6s','6p','5d']
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
      ['Si','Ne',{'3s':2,'3p':2}],\
      ['P' ,'Ne',{'3s':2,'3p':3}],\
      ['S', 'Ne',{'3s':2,'3p':4}],\
      ['Cl','Ne',{'3s':2,'3p':5}],\
      ['Ar','Ne',{'3s':2,'3p':6}],\
      # fourth row
      ['K', 'Ar',{'4s':1,'4p':0}],\
      ['Ti','Ar',{'3d':2,'4s':2,'4p':0}],\
      ['Kr','Ar',{'3d':10,'4s':2,'4p':6}],\
      # fifth row
      ['Sr','Kr',{'5s':2,'4d':0,'5p':0}],
      ['Mo','Kr',{'4d':5,'5s':1,'5p':0}],
      ['Sn','Kr',{'4d':10,'5s':2,'5p':2}],
      ['Pd','Kr',{'4d':10,'5s':0,'5p':0}],
      ['Xe','Kr',{'4d':10,'5s':2,'5p':6}],
      # sixth row
      ['Pt','Xe',{'4f':14,'5d':9,'6s':1,'6p':0}],
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
        print('X'*40)
        for d in data[symbol]:
            print(d,data[symbol][d])

from util import compare_tables
from util import plot_table
from hotbit.parametrization import KSAllElectron
from hotbit.parametrization import SlaterKosterTable

import os

param=os.environ.get('HOTBIT_PARAMETERS')

#Au-Au

#e1=KSAllElectron('Au',confinement={'mode':'quadratic','r0':5.04})
#e1.run()
#e2=e1
#sk=SlaterKosterTable(e1,e2)
#sk.run(1,15,50)
#sk.write()
#compare_tables('Au_Au.par','Au_Au_NR.par',s1='Au',s2='Au',screen=False)
              
lst=[('C','C',1.85*1.46,1.85*1.46),\
     ('C','H',1.85*1.46,1.85*0.705),\
     ('Na','C',1.85*2.9,1.85*1.46),\
     ('O','H',1.85*1.38,1.85*0.705),\
     ('Mg','O',1.85*1.41/0.529177,1.85*1.38),\
     ('Na','O',1.85*2.9,1.85*1.38),\
     ('H','H',1.85*0.705,1.85*0.705)]
     
     
     
for s1,s2,r01,r02 in lst:     
    e1=KSAllElectron(s1,nodegpts=500,confinement={'mode':'quadratic','r0':r01})
    e1.run()
    if s1==s2:  
        e2=e1
    else:   
        e2=KSAllElectron(s2,confinement={'mode':'quadratic','r0':r02})
        e2.run()    
        
    sk=SlaterKosterTable(e1,e2)
    sk.run(1E-3,12,10) #,ntheta=20,nr=20)
    sk.write()
    file='%s_%s.par' %(s1,s2)
    #compare_tables( param+'/'+file,file,s1,s2,screen=False)
    plot_table(file,s1=s1,s2=s2,screen=True,der=0)
    plot_table(file,s1=s1,s2=s2,screen=True,der=1)
        

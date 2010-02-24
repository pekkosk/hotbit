from hotbit.parametrization import SlaterKosterTable
from hotbit.parametrization import KSAllElectron
from util import plot_table
from hotbit import Element
from pickle import load
import pylab as pl
from box.data import data

e1=KSAllElectron('C',confinement={'mode':'quadratic','r0':5.04})
e1.run()
e2=KSAllElectron('H',confinement={'mode':'quadratic','r0':5.04})
e2.run()


sk=SlaterKosterTable(e1,e2)
sk.run(1,15,20,ntheta=50,nr=25)
sk.write()
plot_table('C_H.par',screen=True)   

    
        
        

from hotbit.parametrization import SlaterKosterTable
from hotbit import Element
from pickle import load
#import cgitb; cgitb.enable()
import pylab as pl

e1=Element('H')
e1.read_functions('H.fun')
e1.data['epsilon']={'1s':-0.231746022795451}
e2=e1

sk=SlaterKosterTable(e1,e2)
sk.set_comment('slako comment')
sk.calculate_tables(0.001,3,4)
sk.plot()

    
        
        

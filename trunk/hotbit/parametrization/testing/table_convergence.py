import numpy as nu
import pylab as pl
from hotbit.parametrization import KSAllElectron
from hotbit.parametrization import SlaterKosterTable
s1,s2,r01,r02=('Mg','O',1.85*1.41/0.529177,1.85*1.38)
        
     
     
e1=KSAllElectron(s1,nodegpts=500,confinement={'mode':'quadratic','r0':r01})
e1.run()
if s1==s2:  
    e2=e1
else:   
    e2=KSAllElectron(s2,nodegpts=500,confinement={'mode':'quadratic','r0':r02})
    e2.run()    
    
# searching optimum ntheta/nr for fixed number of grid points    
if False:    
    y=[]
    x=[]
    t=[]
    tlist=[10,20,40,60,100,200,300,400,500,501]
    #rlist=[10,20,40,50,100,200]
    
    
    
    for ntheta in tlist:
        if ntheta!=501:
            nr=10000/ntheta
        else:
            nr=200
        #for nr in rlist:
        #x.append(ntheta*nr)
        #t.append(ntheta)
        sk=SlaterKosterTable(e1,e2)
        sk.run(4,4,1,ntheta=ntheta,nr=nr)
        r,table=sk.get_table()
        y.append([table[0][0,5],table[0][0,6],table[0][0,8],table[1][0,8],table[0][0,9]])
            
    y=nu.array(y)     
    #t=nu.array(t)
    for i in range(5): 
        #pl.scatter(x,y[:,i]-y[-1,i],s=t)        
        pl.plot(tlist,abs(y[:,i]-y[-1,i]))        
    pl.axhline()    
    pl.show()    

# convergence as a number of grid points with fixed ntheta/nr
if True:    
    y=[]
    tlist=[10,20,40,60,100,200,300,400,500]
    for ntheta in tlist:
        nr=ntheta/3
        sk=SlaterKosterTable(e1,e2)
        sk.run(4,4,1,ntheta=ntheta,nr=nr)
        r,table=sk.get_table()
        y.append([table[0][0,5],table[0][0,6],table[0][0,8],table[1][0,8],table[0][0,9]])
            
    y=nu.array(y)     
    for i in range(5): 
        pl.semilogy(tlist,abs(y[:,i]-y[-1,i])+1E-9)        
    pl.axhline()    
    pl.show()    

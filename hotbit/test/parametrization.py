from hotbit.parametrization import KSAllElectron, SlaterKosterTable

atom=KSAllElectron('K', xc='PW92',
                    confinement={'mode':'quadratic','r0':12.})
#                   confinement={'mode':'Woods-Saxon','r0':8.,'a':4.,'W':0.75})
atom.run()
atom.plot_Rnl()

table=SlaterKosterTable(atom,atom)
table.run(R1=1,R2=10,N=3,ntheta=50,nr=10)        
table.plot()    
        
    
        

from hotbit.parametrization import KSAllElectron, SlaterKosterTable

atom=KSAllElectron('C',xc='pw92')#,confinement={'mode':'quadratic','r0':5.3})#,txt='-')
atom.run()
atom.plot_Rnl()

#table=SlaterKosterTable(atom,atom)#,txt='-')
#table.run(R1=1,R2=10,N=3,ntheta=50,nr=10)        
#table.plot()    
        
    
        

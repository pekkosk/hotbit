from hotbit.parametrization import KSAllElectron, SlaterKosterTable

atom=KSAllElectron('C',convergence={'density':1E-1,'energies':1E-1},txt='/dev/null')
atom.run()
    
table=SlaterKosterTable(atom,atom,txt='/dev/null')
table.run(R1=1,R2=10,N=3,ntheta=50,nr=10,wflimit=1E-7)        
    
        
    
        

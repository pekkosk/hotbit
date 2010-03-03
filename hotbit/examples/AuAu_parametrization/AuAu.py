from hotbit import * 
from ase import *

phase = 3

if phase==1:
    # solve free atoms 
    # Make Au.elm:
    # * on-site -energies
    # * calculate or find values for U, FWHM
    Au = KSAllElectron('Au',scalarrel=True,txt='Au_free.txt')    
    Au.run()
    Au.plot_Rnl('Au_free.pdf')
    
    
    
if phase==2:
    # solve confined atoms
    r0 = 1.85 * 1.36/Bohr
    Au = KSAllElectron('Au',txt='Au.txt',confinement={'mode':'quadratic','r0':r0})    
    Au.run()
    Au.plot_Rnl('Au.pdf')
    

    # calculate Slater-Koster tables and output
    table = SlaterKosterTable(Au,Au)
    table.run(1,12,50)
    table.write('Au_Au_no_repulsion.par')
    table.plot('Au_Au_slako.pdf')


if phase==3:
    # fit repulsion
    tab={'AuAu':'Au_Au_no_repulsion.par'}
    elm={'Au':'Au.elm'}
    mixer={'name':'Anderson','mixing_constant':0.1,'convergence':1E-9}
    calc0 = Hotbit(txt='-',elements=elm,mixer=mixer,tables=tab,SCC=True)
    calc1 = Hotbit(txt='-',elements=elm,mixer=mixer,tables=tab,SCC=False,kpts=(6,6,6))
    calc2 = Hotbit(txt='-',elements=elm,mixer=mixer,tables=tab,SCC=True,charge=-1.0)
    
    rep = RepulsiveFitting('Au','Au',r_cut=3.3,s=100)
    
    # add data
    #rep.append_point(weight=0.1,R=3.2,dvrep=-0.2,comment='example point by hand')
    rep.append_dimer(weight=0.5,calc=calc0,R=2.49,comment='Au2')
    rep.append_energy_curve(weight=1.0,calc=calc0,traj='dimer_curve.traj',label='DFT dimer',comment='dimer curve')
    rep.append_scalable_system(weight=1.0,calc=calc1,atoms='bulk.traj',comment='bulk Au',label='bulk')
    rep.append_homogeneous_cluster(weight=1.0,calc=calc2,atoms='Au12-.xyz')
    
 
    # fit & output
    rep.fit()
    rep.add_comment('Repulsion by Pekka Koskinen')
    rep.write_par('Au_Au_no_repulsion.par',filename='Au_Au_repulsion.par')
    rep.plot('AuAu_repulsion.pdf')
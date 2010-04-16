from numpy import *
from hotbit import * 
from ase import *

phase = 3

if phase==1:
    # solve free atoms 
    # Make C.elm:
    # * on-site -energies
    # * calculate or find values for U, FWHM
    A = KSAllElectron('C',txt='C_free.txt')    
    A.run()
    A.plot_Rnl('C_free.pdf')
    
    B = KSAllElectron('H',txt='H_free.txt')    
    B.run()
    B.plot_Rnl('H_free.pdf')
    
    
    
if phase==2:
    # solve confined atoms
    r0A = 1.85 * 0.76/Bohr
    A = KSAllElectron('C',txt='C.txt',confinement={'mode':'quadratic','r0':r0A})    
    A.run()
    A.plot_Rnl('C.pdf')
    A.write_unl('C.wf')
    
    r0B = 1.85 * 0.31/Bohr
    B = KSAllElectron('H',txt='H.txt',confinement={'mode':'quadratic','r0':r0B})    
    B.run()
    B.plot_Rnl('H.pdf')
    B.write_unl('H.wf')
    
    # calculate Slater-Koster tables and output
    table = SlaterKosterTable(A,B)
    min, rmax, dr = 1,10,0.1
    N = int( (rmax-rmin)/dr )
    print N,'points in the table'
    table.run(rmin,rmax,N)
    table.write('C_H_norep.par')
    table.plot('C_H_slako.pdf')


if phase==3:
    # fit repulsion
    tab={'CH':'C_H_norep.par','others':'default'}
    elm={'C':'C.elm','H':'H.elm'}
    mixer={'name':'Anderson','mixing_constant':0.1,'convergence':1E-9}
    calc0 = Hotbit(width=0.1,txt='-',elements=elm,tables=tab,mixer=mixer,SCC=True)
    calc1 = Hotbit(width=0.1,txt='-',elements=elm,tables=tab,mixer=mixer,SCC=True,charge=-1)
    
    rep = RepulsiveFitting('C','H',r_cut=1.6,s=1.0)
    
    # add data
    rep.append_dimer(weight=1.0,calc=calc0,R=1.137,comment='CH dimer')
    rep.append_energy_curve(weight=1.0,calc=calc1,traj='CH-.traj',label='CH-',comment='CH-')
    rep.append_energy_curve(weight=1.0,calc=calc0,traj='ethyne.traj',comment='ethyne')
    rep.append_energy_curve(weight=1.0,calc=calc0,traj='methane.traj',comment='methane')
    rep.append_energy_curve(weight=1.0,calc=calc0,traj='benzene.traj',comment='benzene')
    rep.write_fitting_data('fitting_data.dat')       
 
    # fit & output
    rep.fit()
    rep.add_comment('Repulsion by Pekka Koskinen')
    rep.write_par('C_H_norep.par',filename='C_H_rep.par')
    rep.plot('CH_rep.pdf')
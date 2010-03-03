from hotbit import * 
from ase import *

phase = 3

if phase==1:
    # solve free atoms 
    # Make C.elm:
    # * on-site -energies
    # * calculate or find values for U, FWHM
    C = KSAllElectron('C',txt='C_free.txt')    
    C.run()
    C.plot_Rnl('C_free.pdf')
    
    H = KSAllElectron('H',txt='H_free.txt')    
    H.run()
    H.plot_Rnl('H_free.pdf')
    
    
    
if phase==2:
    # solve confined atoms
    r0C = 1.85 * 0.76/Bohr
    C = KSAllElectron('C',txt='C.txt',confinement={'mode':'quadratic','r0':r0C})    
    C.run()
    C.plot_Rnl('C.pdf')
    
    r0H = 1.85 * 0.31/Bohr
    H = KSAllElectron('H',txt='H.txt',confinement={'mode':'quadratic','r0':r0H})    
    H.run()
    H.plot_Rnl('H.pdf')
    
    # calculate Slater-Koster tables and output
    table = SlaterKosterTable(C,H)
    table.run(1,12,50)
    table.write('C_H_no_repulsion.par')
    table.plot('C_H_slako.pdf')


if phase==3:
    # fit repulsion
    tab={'CH':'C_H_no_repulsion.par','others':'default'}
    elm={'C':'C.elm','H':'H.elm'}
    mixer={'name':'Anderson','mixing_constant':0.1,'convergence':1E-9}
    calc0 = Hotbit(width=0.1,txt='-',elements=elm,tables=tab,mixer=mixer,SCC=True)
    calc1 = Hotbit(width=0.1,txt='-',elements=elm,tables=tab,mixer=mixer,SCC=True,charge=-1)
    
    rep = RepulsiveFitting('C','H',r_cut=1.6,s=100)
    
    # add data
    rep.append_dimer(weight=1.0,calc=calc0,R=1.137,comment='CH dimer')
    rep.append_energy_curve(weight=1.0,calc=calc1,traj='CH-.traj',label='CH-',comment='CH-')
    rep.append_energy_curve(weight=1.0,calc=calc0,traj='ethyne.traj',label='ethyne',comment='ethyne')
    rep.append_energy_curve(weight=1.0,calc=calc0,traj='methane.traj',label='methane',comment='methane')
    rep.append_energy_curve(weight=1.0,calc=calc0,traj='benzene.traj',label='benzene',comment='benzene')
    rep.write_fitting_data('fitting_data.dat')       
 
    # fit & output
    rep.fit()
    rep.add_comment('Repulsion by Pekka Koskinen')
    rep.write_par('C_H_no_repulsion.par',filename='C_H_repulsion.par')
    rep.plot('CH_repulsion.pdf')
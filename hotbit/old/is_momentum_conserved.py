#!/usr/bin/env python
"""
    Inspect the conservation of total momentum in MD run.
    
    Read electron and ion momenta from "momenta.out", calculate
    their variances and compare.
    
    P. Koskinen 29.3 2007
    
"""
import box.mix as mix
import numpy,os

dat = mix.read('momenta.out','r')
Pix = dat[:,0]
Piy = dat[:,1]
Piz = dat[:,2]
Pex = dat[:,3]
Pey = dat[:,4]
Pez = dat[:,5]
Px  = Pix+Pex
Py  = Piy+Pey
Pz  = Piz+Pez
av  = numpy.average
va  = numpy.var
print '-------------------------------------------------'
print 'Conservation of momentum:'
print 'Averages for P (ions) =', av(Pix),av(Piy),av(Piz)
print 'Averages for P (elec) =', av(Pex),av(Pey),av(Pez)
print 'Averages for P (tot)  =', av(Px),av(Py),av(Pz)
print 'Variance for P (ions) =', va(Pix),va(Piy),va(Piz)
print 'Variance for P (elec) =', va(Pex),va(Pey),va(Pez)
print 'Variance for P (tot)  =', va(Px),va(Py),va(Pz)
print '-------------------------------------------------'

plot="""set term pos col enhanced "Helvetica" 18
    set output 'conserv_P.ps'
    set title 'Testing dftb; check conservation of momentum'
    set xlabel 'time step'
    set ylabel 'P_i (a.u.)'
    pl 'momenta.out' u ($1+$4) w l t 'P_x',\
       '' u ($2+$5) w l t 'P_y',\
       '' u ($3+$5) w l t 'P_z'"""
mix.gplot(plot,'conserv_P.gp')
psviewer = os.environ.get('PSVIEWER')
mix.execute('%s conserv_P.ps' %psviewer )
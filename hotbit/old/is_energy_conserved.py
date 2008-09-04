#!/usr/bin/env python
"""
    Inspect the conservation of total energy in microcanonical run.
    
    Read total energy ('etot') and potential energy ('epot')
    from the file 'loop.out' and calculate their variances
    and compare.
    
"""
import box.mix as mix
import numpy as N
import os

dat = mix.read('loop.out','r')
f   = open('loop.out','r')
list= f.readline().lstrip('#').split()
f.close()
for i in range(len(list)):
    if( list[i].lower()=='etot' ): ti=i
    if( list[i].lower()=='epot' ): pi=i
    
etot = N.array(dat[:,ti])
epot = N.array(dat[:,pi])

print '-------------------------------------------------'
print 'Conservation of energy:'
dE=N.sqrt(N.var(etot))
dU=N.sqrt(N.var(epot))
print 'deviation for total energy    =', dE,'(',dE*27.2114*1000,'meV)'
print 'deviation for potential energy=', dU
print 'sigma(etot)/sigma(epot)      =',dE/dU
print '-------------------------------------------------'

pi+=1
ti+=1
plot="""set term pos col enhanced
    set output 'const_E.ps'
    set title 'testing dftb; constant-energy run'
    pl 'loop.out' u 1:%i w l t 'E_{pot}','' u 1:%i w l t 'E_{tot}' """ %(pi,ti)
mix.gplot(plot,'const_E.gp')
psviewer = os.environ.get('PSVIEWER')
mix.execute('%s const_E.ps' %psviewer )
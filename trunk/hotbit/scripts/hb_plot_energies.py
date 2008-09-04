#!/usr/bin/env python
import box.mix as mix
import os
exists=mix.find_value('md.dat','out_loop',fmt='test')
if exists: line=mix.find_value('md.dat','out_loop','all')
else:      line=mix.find_value('source/defaults.dat','out_loop','all')

l={'etot':1,'epot':2}
for i in range(len(line)):
    item = line[i]
    if l.has_key(item): l[item]=i

etot=mix.read_column('etot','loop.out')
e_0 =etot[0]
emax=max(etot)
emin=min(etot)

eti =l['etot']+1
epi =l['epot']+1
file='plot_energies'
plot="""set term pos color enhanced 'Helvetica' 20
    set output '%(file)s.ps'
    E_0=%(e_0)f
    HA=27.2114
    #set yra [:0.1]
    set ylabel 'energy (eV)'
    set xlabel 'output step'
    set title 'Total, potential etc. energy info'
    pl 'loop.out' u 0:(($%(eti)i-E_0)*HA) w l t 'etot',\
               '' u 0:(($%(epi)i-E_0)*HA) w l t 'epot'""" %vars()

print 'max-min for etot',(emax-emin)*27.2114,'eV'

mix.gplot(plot,file+'.gp')
psviewer = os.environ.get('PSVIEWER')
mix.execute('%s %s.ps &' %(psviewer,file) )

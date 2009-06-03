#!/usr/bin/env python
"""
    Plot the energies and occupations.
    
    Usage:
        hb_DOS.py file
    
    Example:
        hb_DOS.py dos
        
    P. Koskinen 19.4 2007
"""
import box.mix as mix
import sys,os
import numpy as N

try:
    file=sys.argv[1]
except:
    print __doc__
    sys.exit(1)    

e=mix.read('ev.out')
ef=max( N.compress(N.greater(e[:,2],1E-4),e[:,1]) )
plot="""set term pos color enhanced 'Helvetica' 24
    set output '%s.ps'
    HA=27.2114
    EF=%f
    set title 'States and occupations'
    set xlabel 'energy (eV)'
    set ylabel 'occupations'
    set yrange [:4]
    pl 'ev.out' u (($2-EF)*HA):($3) w i t 'occupied',\\
             '' u (($2-EF)*HA):(2-$3) w i t 'unoccupied'""" %(file,ef)
    
mix.gplot(plot,file+'.gp')
psviewer=os.environ.get('PSVIEWER')
mix.execute('%s %s.ps' %(psviewer,file) )



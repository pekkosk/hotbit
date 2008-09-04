#!/usr/bin/env python
"""
    Make gnuplot file and plot the linear optical response from the file optical.out.
    
    Usage:
        hb_optical.py emax
        
        emax - the maximum energy in eV plotted
        
    Example:
        hb_optical.py 5.0
        
    P. Koskinen 21.3 2007
"""
import box.mix as mix
import os,sys

try:
    emax=float(sys.argv[1])
except:
    print __doc__
    sys.exit(1)

plot="""set term pos color enhanced "Helvetica" 20
        set output 'optical.ps'
        set key off
        set title 'Linear optical response'
        set xlabel 'energy (eV)'
        set ylabel 'Lorenzian broadened response'
        set xra [0:%(emax)f]
        pl 'optical.out' i 0 u ($1*27.2114):2 w l """ %vars()
mix.gplot(plot,'optical.gp')
psviewer=os.environ.get('PSVIEWER')
mix.execute('%s optical.ps' %psviewer)
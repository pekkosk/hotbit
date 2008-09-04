#!/usr/bin/env python
"""
    Plot energy occupations for monitoring level structures and hoppings.
    
    Pekka Koskinen 24.5 2007
"""
from box import mix
import os

f=open('ev.out','r')
traj=[]
while 1:
    frame=mix.read(f)
    if frame==None: break
    traj.append(frame)
f.close()

F=len(traj)
o=open('plot_occu.out','w')
# data for showing occupations in energy ordering
for i in range(F):
    frame=traj[i]
    for (k,f) in zip(frame[:,0],frame[:,2]):
        print>>o,i,k,f
        
o.write('\n\n')
# data for showing occupations in energy scale
for i in range(F):
    frame=traj[i]
    for (e,f) in zip(frame[:,1],frame[:,2]):
        print>>o,i,e,f
o.close()
        
       
plot="""set term pos color enhanced "Helvetica" 10
        set output 'plot_occu.ps'
        
        set size 1,1
        set multiplot
        HA=27.2114
        set pointsize 0
        
        set size 0.5,1.0
        set origin 0,0
        set key off
        set xlabel 'MD output step'
        set ylabel 'occupations for energy-ordered levels'
        pl 'plot_occu.out' i 0 u 1:($2) w p ps 1 pt 1,\
                        '' i 0 u 1:(($3)>=0.99 ? $2:1/0) w p ps 1 pt 2,\
                        '' i 0 u 1:(($3)>=1.99 ? $2:1/0) w p ps 2 pt 6
                        
                        
        
        set size 0.5,1.0
        set origin 0.5,0
        set key off
        set xlabel 'MD output step'
        set ylabel 'occupations on energy scale (eV)'
        pl 'plot_occu.out' i 1 u 1:(($2)*HA) w p ps 1 pt 1,\
                        '' i 1 u 1:(($3)>=0.99 ? ($2)*HA:1/0) w p ps 1 pt 2,\
                        '' i 1 u 1:(($3)>=1.99 ? ($2)*HA:1/0) w p ps 2 pt 6
    """
        

mix.gplot(plot,'plot_occu.gp')
viewer=os.environ.get('PSVIEWER')
mix.execute('%s plot_occu.ps &' %viewer)
        
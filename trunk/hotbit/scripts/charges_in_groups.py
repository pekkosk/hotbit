#!/usr/bin/env python
"""
    Calculate total charges in different atom groups.
    
    Look at the beginning of the script for default input and
    output files.
    
    Usage:
        charges_in_groups.py         
    
    P. Koskinen 27.3 2007
"""
import box.mix as mix
import box.md as md
import numpy as nu

atoms='atoms.dat'
charges='dq.out'
out='dq_groups.out'


mol=md.Molecule(atoms,format='dat')
f=open(charges,'r')
N=mol.N
dq_grp= [[] for i in range(N)]
ind   = [[] for i in range(N)]
o=open(out,'w')
while 1:
    dq=mix.read(f)
    if dq==None or len(dq[:])<N: break
    dq_g= nu.zeros(N)
    sum = [0 for i in range(N)]
    Ng  = 0
    for i in range(N):
        g = mol.atoms[i].act-1
        if sum[g]==0: Ng=Ng+1
        sum[g]+=1
        dq_g[g]+=float(dq[i])
    for g in range(N): dq_grp[g].append(dq_g[g])

print "Number of groups",Ng        
for g in range(N):
    if sum[g]==0: continue
    o.write('#group %i\n' %(g+1))
    for dq in dq_grp[g][:]: o.write('%f\n' %dq)
    o.write('\n\n')
       
f.close()

#mix.gplot(plot,'dq_groups.gp')
#viewer=os.penviron.get('PSVIEWER')
#os.system('%s dq_groups.ps' %viewer)

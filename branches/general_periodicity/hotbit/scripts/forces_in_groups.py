#!/usr/bin/env python
"""
    Calculate total forces acting in different atom groups.
    
    Look at the beginning of the script for default input and
    output files.
    
    Usage:
        forces_in_groups.py         
    
    P. Koskinen 30.3 2007
"""
import box.mix as mix
import box.md as md
import numpy as nu

atoms  ='atoms.dat'
forces ='forces.out'
out    ='forces_groups.out'


mol=md.Molecule(atoms,format='dat')
file=open(forces,'r')
N=mol.N
F_grp=[[] for i in range(N)] #forces for all groups and time steps
o=open(out,'w')
while 1:
    F=mix.read(file)
    if len(F)<N: break
    flag = [0 for i in range(N)]
    F_g  = [nu.array([0.,0.,0.]) for i in range(N)]
    for i in range(N):
        g = mol.atoms[i].act
        flag[g]+=1
        F_g[g][:]+=F[i,:]
    for g in range(N): F_grp[g].append(F_g[g][:])
        
for g in range(N):
    if flag[g]==0: continue
    o.write('#forces in group %i\n' %g)
    for F in F_grp[g]: o.write('%f %f %f\n' %(F[0],F[1],F[2]))
    o.write('\n\n')

file.close()

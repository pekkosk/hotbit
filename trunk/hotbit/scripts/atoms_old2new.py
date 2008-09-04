#!/usr/bin/env python
"""
    Transform atoms.dat in the old format ("fixed format"-like)
    into the new format (free, keywork-based format)
    
    Example:
        atoms_old2new.py atoms_old.dat atoms_new.dat
        
    P. Koskinen 6.3 2007
"""
import os,sys

try:
    inp  = sys.argv[1]
    outp = sys.argv[2]
except:
    print __doc__
    sys.exit(1)


#
# read the old format
#
fi=open(inp,'r')
fi.readline()
N=int(fi.readline())
fi.readline()
occ=float(fi.readline())
fi.readline()
atoms=[]; v=[]; f=[]
#el,m,x,y,z,act,gamma,Tlang
for i in range(N): atoms.append( fi.readline().split() )
# velocities    
fi.readline()
for i in range(N): v.append( fi.readline() )
# forces    
fi.readline()
for i in range(N): f.append( fi.readline() )
fi.close()

#
# write in the new format
#
fo=open(outp,'w')
fo.write('nat=%i\n' %N)
fo.write('extra_electrons=[occ=%f in the old format]\n' %occ)
fo.write('atoms=\n')
for i in range(N): 
    el,x,y,z=(atoms[i][0],atoms[i][2],atoms[i][3],atoms[i][4])
    #x=str(float(x)/0.529177) # for files with Ångström
    #y=str(float(y)/0.529177)
    #z=str(float(z)/0.529177)
    fo.write( '%s %s %s %s\n' %(el,x,y,z) )
fo.write('flags=\n')
for i in range(N): fo.write(atoms[i][5]+'\n')
fo.write('velocities=\n')
for vel in v: fo.write(vel)
fo.write('forces=\n')
for force in f: fo.write(force)
fo.write('langevin=\n')
for i in range(N): fo.write(atoms[i][6]+' '+atoms[i][7]+'\n')
fo.close()


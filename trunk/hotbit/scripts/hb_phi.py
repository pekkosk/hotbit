#!/usr/bin/env python
"""
    Output the (internal) electrostatic potential into a grid data for vtk file.
    
    Usage:
        hb_phi.py
        phi_tb.out --> phi_tb.vtk
    
    P. Koskinen 126.4 2007
"""
import box.mix as mix 

fi=open('phi_tb.out')
grids=mix.find_value(fi,'xyz_grid',fmt='strings')
data=mix.find_value(fi,'data',fmt='strings')

nx=len(grids[0].split())
ny=len(grids[1].split())
nz=len(grids[2].split())

of=open('phi_tb.vtk','w')
print>>of,"# vtk DataFile Version 2.0"
print>>of,"NOTB simulation"
print>>of,"ASCII"
print>>of,"DATASET RECTILINEAR_GRID"
print>>of,"DIMENSIONS %i %i %i" %(nx,ny,nz)
print>>of,"X_COORDINATES %i double" %nx
print>>of,grids[0]
print>>of,"Y_COORDINATES %i double" %ny
print>>of,grids[1]
print>>of,"Z_COORDINATES %i double" %nz
print>>of,grids[2]
print>>of,"POINT_DATA %i" %(nx*ny*nz)
print>>of,"SCALARS phi double"
print>>of,"LOOKUP_TABLE default"
mx=-1000
mn=1000
for line in data:
    mn=min(mn,float(line))
    mx=max(mx,float(line))
    print>>of,line

print 'min ... max=',mn,'...',mx
of.close()
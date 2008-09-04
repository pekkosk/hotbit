#!/usr/bin/env python
"""
    Output wave functions into a grid data for vtk file.
    
    Transforms files plot_wf_X.out info wf_X.vtk.
    Looks for the state indices ('X') from md.dat (=plot_wf_list)
       
    Usage:
        hb_wf.py
    
    P. Koskinen 27.6 2007
"""
import box.mix as mix 

print 'Make vtk files for plotting for the selected wave functions'
lst=mix.find_value('md.dat','plot_wf_list',fmt='all')

for state in lst:
    print 'plot state wf_%s.out' %state
    fi=open('wf_%s.out' %state)
    msteps=1000
    for step in range(1,msteps+1):
        exist=mix.find_value(fi,'xyz_grid_%i' %step,fmt='test')
        if not exist: break
        grids=mix.find_value(fi,'xyz_grid_%i' %step,fmt='strings')
        data=mix.find_value(fi,'data_%i' %step,fmt='strings')
        
        nx=len(grids[0].split())
        ny=len(grids[1].split())
        nz=len(grids[2].split())
        
        of=open('wf_%s_%i.vtk' %(state,step-1),'w')
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
        print>>of,"SCALARS wavef double"
        print>>of,"LOOKUP_TABLE default"
        mx=-1000
        mn=1000
        for line in data:
            mn=min(mn,float(line))
            mx=max(mx,float(line))
            print>>of,line
        
        print 'min ... max=',mn,'...',mx
        of.close()
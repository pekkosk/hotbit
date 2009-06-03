import box.mix as mix 

def rectilinear_vtk(grid,data,fname):
    """ Write data in rectilinear grid into .vtk file. 
    
    parameters:
    -----------
    grid: grid[:,{0,1,2}] x-, y-, and z-grids
    data: data on this regular grid.
    fname: output file name
    """
    f=open(fname,'w')    
    nx, ny, nz=len(grid[0][:]), len(grid[1][:]), len(grid[2][:])
    print>>f,"# vtk DataFile Version 2.0"
    print>>f,"...some rectilinear grid data."
    print>>f,"ASCII"
    print>>f,"DATASET RECTILINEAR_GRID"
    print>>f,"DIMENSIONS %i %i %i" %(nx,ny,nz)
    print>>f,"X_COORDINATES %i double" %nx
    print>>f,mix.a2s(grid[0][:])
    print>>f,"Y_COORDINATES %i double" %ny
    print>>f,mix.a2s(grid[1][:])
    print>>f,"Z_COORDINATES %i double" %nz
    print>>f,mix.a2s(grid[2][:])
    print>>f,"POINT_DATA %i" %(nx*ny*nz)
    print>>f,"SCALARS some_data double"
    print>>f,"LOOKUP_TABLE default"
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                print>>f, data[i,j,k]
    print 'min ... max=',min(data.flatten()),'...',max(data.flatten())
    f.close()

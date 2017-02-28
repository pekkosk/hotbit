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



def atoms_vtk(atoms,scalars={},vectors={},filename=None):
    '''
    vtk output of atoms
         
    @param filename: vtk output file name
    @param atoms:    atoms object
    @param scalars:  dictionary of atoms' scalar properties
    @param vectors:  dictionary of atoms' vector properties
    '''
    if filename is None:
        filename=atoms.get_name()+'.vtk'
    N=len(atoms)
    f=open(filename,'w')
    f.write('# vtk DataFile Version 2.0 \nAtoms %s\n' %atoms.get_name())
    f.write('ASCII \nDATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS %i double \n ' %N)
    fmt='%20.14f' #output format for floats
    
    # Point data (atom coordinates) and cell data (bonds)
    for r in atoms.get_positions():
        f.write('%s\n' %mix.a2s(r,fmt=fmt))
        
    # First the data related to atoms
    f.write('POINT_DATA %i\n' %N)
    for scalar in scalars:
        print>>f, 'SCALARS %s double 1\nLOOKUP_TABLE default' %scalar
        for value in scalars[scalar]:
            print>>f, '%12.6f' %(value*1.0)
    for vector in vectors:
        print>>f, 'VECTORS %s double\n' %vector
        for value in properties:
            print>>f, mix.a2s(value,fmt=fmt)
    f.close()     

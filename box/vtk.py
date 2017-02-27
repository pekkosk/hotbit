from __future__ import print_function

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
    print("# vtk DataFile Version 2.0", file=f)
    print("...some rectilinear grid data.", file=f)
    print("ASCII", file=f)
    print("DATASET RECTILINEAR_GRID", file=f)
    print("DIMENSIONS %i %i %i" %(nx,ny,nz), file=f)
    print("X_COORDINATES %i double" %nx, file=f)
    print(mix.a2s(grid[0][:]), file=f)
    print("Y_COORDINATES %i double" %ny, file=f)
    print(mix.a2s(grid[1][:]), file=f)
    print("Z_COORDINATES %i double" %nz, file=f)
    print(mix.a2s(grid[2][:]), file=f)
    print("POINT_DATA %i" %(nx*ny*nz), file=f)
    print("SCALARS some_data double", file=f)
    print("LOOKUP_TABLE default", file=f)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                print(data[i,j,k], file=f)
    print('min ... max=',min(data.flatten()),'...',max(data.flatten()))
    f.close()



def atoms_vtk(atoms,scalars={},vectors={},filename=None):
    '''
    vtk output of atoms
         
    @param filename: vtk output file name
    @param atoms:    atoms object
    @param scalars:  dictionary of atoms' scalar properties
    @param vectors:  dictionary of atoms' vector properties
    '''
    if filename==None:
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
        print('SCALARS %s double 1\nLOOKUP_TABLE default' %scalar, file=f)
        for value in scalars[scalar]:
            print('%12.6f' %(value*1.0), file=f)
    for vector in vectors:
        print('VECTORS %s double\n' %vector, file=f)
        for value in properties:
            print(mix.a2s(value,fmt=fmt), file=f)
    f.close()     
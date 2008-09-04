#!/usr/bin/env python
"""
    Make VTK folder and data from the available HOTBIT simulation data.
    
    The output will be in ./vtk-folder
    
    Usage:
        
        hb_vtk.py S ["add1 x y z" "add2 x y z"]
        
        * 'S' skip every S'th md step
        
        * '"add x y z"' *(note: in quotation marks)*
        
            Data can be supplemented by additional scalar or vector data. 
            Available options are: 
                
                * 'ekin'
                
                * 'epot'
                
                * 'dipole'
                
                * 'E_field'
                
            which will be shown in given points (x,y,z) in space.
           
    examples:        
        
        * hb_vtk.py 1 
        
            makes vtk files for viewing trajectory, forces, charges etc.
            
        * hb_vtk.py 1 "dipole 10 0 0"
        
            as above, but appends the dipole moment arrow at (10,0,0)
    
    Pekka Koskinen 27.4 2007
"""

import box.mix as mix
import box.hotbit_tools as hotbit_tools
import box.md as md
import numpy as nu
from os.path import isfile as opi
import sys,os

file_atoms      ='atoms.dat'
file_traj       ='traj.xyz'
file_charges    ='dq.out'
file_forces     ='forces.out'
file_f_sep      =['forces_1.out','forces_2.out','forces_3.out']
file_velocities ='velocities.out'
file_bonds      ='bonding.out'
file_momenta    ='momenta.out'
file_loop       ='loop.out'
file_dipole     ='dipolemom.out'
file_E_field    ='E_field.out'
file_vectors    ='vectors.out';  vectors_name='NA_coupling'
vector_files    = {'nac_2_1':'nac_2_1.out',
                   'nac_3_1':'nac_3_1.out',
                   'nac_3_2':'nac_3_2.out',
                   'nac_1_2':'nac_1_2.out'}
output          ='vtk_output'

try:
    S = int(sys.argv[1])
except:
    print __doc__
    sys.exit(1)

N_add=len(sys.argv[1:])-1
name_list=[]
points=[]
if os.path.isdir(output): 
    ans=raw_input('%s exists. Remove the whole directory? (yes/no) ' %output)
    if ans=='yes': mix.execute('rm -r %s' %output)
    else: 
        print 'Abort.'
        sys.exit(1)
    
if N_add!=0:
    toadd=sys.argv[2:]
    for item in toadd:
        try:    name,x,y,z = item.split()
        except: mix.error_exit('Invalid input arguments.')
        print 'Adding',name,'to position',x,y,z
        name_list.append( name )
        points.append( [float(x),float(y),float(z)] )
        
   
#--------------------------------------------
# Main vtk trajectory with atom-specific data
#--------------------------------------------
traj = hotbit_tools.read_trajectory(file_traj)

if opi(file_atoms):         hotbit_tools.read_flags(traj,file_atoms)
if opi(file_charges):       hotbit_tools.read_charges(traj,file_charges)
if opi(file_forces):        hotbit_tools.read_forces(traj,file_forces)
if opi(file_velocities):    hotbit_tools.read_velocities(traj,file_velocities)
if opi(file_bonds):         hotbit_tools.read_bonds(traj,file_bonds)
for item in vector_files:
    file=vector_files[item]
    if opi(file):
        hotbit_tools.read_vectors(traj,file,item)
for i in range(len(file_f_sep)):
    file=file_f_sep[i]
    if opi(file): hotbit_tools.read_vectors(traj,file,'forces_%i' %i)
hotbit_tools.traj_vtk_output(traj,'traj',output,S)





#--------------------------------------
# Additional info in separate vtk files
#--------------------------------------
loop=mix.read('loop.out')
tsteps=len(loop[:,0])

# read all the data
epi = mix.identify_column('epot',file_loop)
eki = mix.identify_column('ekin',file_loop)
if opi(file_dipole):  dip = mix.read(file_dipole)
if opi(file_E_field): E   = mix.read(file_E_field)

#--------------------------
# collect the selected data 
#--------------------------
data=[]
for i in range(tsteps):
    if nu.mod(i,S)!=0: continue
    d=[]
    ind=0
    for name in name_list:
        if name=='epot': d.append( ['scalar','Epot',ind,loop[i,epi]] )
        if name=='ekin': d.append( ['scalar','Ekin',ind,loop[i,eki]] )
        if name=='dipole': d.append( ['vector','Dipole_moment',ind,dip[i,1],dip[i,2],dip[i,3]] )
        if name=='E_field': d.append( ['vector','Electric_field',ind,E[i,0],E[i,1],E[i,2]] )
        ind+=1
    data.append(d)
    
#--------------------------
#    output into *.vtk
#--------------------------
for k in range(tsteps):
    i=k/S
    if nu.mod(k,S)!=0: continue
    
    # initialize
    f=open('%s/add_%i.vtk' %(output,i),'w')
    f.write('# vtk DataFile Version 2.0 \nNOTB simulation\n')
    f.write('ASCII \nDATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS %i double \n ' %N_add)
    for p in points:
        f.write('%f %f %f\n' %(p[0],p[1],p[2]) )
    f.write('POINT_DATA %i\n' %N_add)
    
    
    # output of different data
    for d in data[i]:
        if d[0]=='scalar':
            f.write('SCALARS %s double 1\nLOOKUP_TABLE default\n' %d[1])       
            for j in range(N_add):
                if j!=d[2]: f.write('0\n')
                else: f.write('%f\n' %d[3])
        if d[0]=='vector':
            f.write('VECTORS %s double\n' %d[1])
            for j in range(N_add):
                if j!=d[2]: f.write('0 0 0\n')
                else: f.write('%f %f %f\n' %(d[3],d[4],d[5]))
    f.close()

        

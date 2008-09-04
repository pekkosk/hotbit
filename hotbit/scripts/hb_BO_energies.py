#!/usr/bin/env python
"""
    Calculate the Born-Oppenheimer energies for a given MD trajectory.
    
    Useful for calculating the electronic excitation energy in TDTB calculations
    compared to BO surface.
    
    Usage:
        hb_BO_energies.py n traj loop
     
        n - calculate BO energy for every n'th frame
        traj - the trajectory (xyz) file (e.g. from TDTB calculation)
        loop - the loop output file (e.g. from TDTB calculation)
        
    Example:
        hb_BO_energies.py 10 traj_TD.xyz loop_TD.out
        
    Note: 
        - The output for trajectory file and loop-file has to have the same frequency
          (out_freq1 has to equal to out_freq2)
        - The current directory should have a HOTBIT code ready for communication.
        - output is a file energy_BO.out listing first the time step and the BO energies,
          and later the differences of TDTB and BO energies
    
    P. Koskinen 5.4 2007        
"""
import box.mix as mix
import box.md as md
import box.hotbit_tools as hotbit_tools
import numpy as N
import sys,os

try:
    every,ftraj,floop =sys.argv[1:]
except:
    print __doc__
    sys.exit(0)

every   =float(every)
ofile   ='energy_BO.out'
tb_code ='.'

com=mix.find_value('md.dat','communicate')
if com!='yes': mix.error_exit('hotbit code is not communicating!')

traj=hotbit_tools.read_trajectory(ftraj)
#traj[0].output_atoms_dat('atoms.dat')
print len(traj), "MD steps, calculate BO energy for every", every, "steps"
i=-1
of=open(ofile,'w')
print>>of,"#md step, BO energy"

# calculate BO energies for the trajectory
E_BO=[]
steps=[]
mol=md.Molecule('atoms.dat')
files=hotbit_tools.hotbit_interactive(mol,'start')
for mol in traj:
    i+=1
    if N.mod(i,every)!=0: continue
    #epot=hotbit_tools.hotbit_energy(mol,tb_code,'Epot')
    epot=hotbit_tools.hotbit_interactive(mol,'energy',files)
    print>>of, i+1, epot    
    print i+1,epot
    E_BO.append(epot)
    steps.append(i+1)

print>>of
print>>of
ret=hotbit_tools.hotbit_interactive(mol,'stop',files)


# Compare with potential energies 
# from a TD-TB calculation
print>>of,"#md step, TD energy above BO energy"
energies_TD=mix.read_column('epot',floop)
E_TD=[]
i=-1
for e in energies_TD:
    i+=1
    if N.mod(i,every)!=0: continue
    E_TD.append(e)
    
if len(E_TD)!=len(E_BO): 
    print 'BO and TD have different lengths:'
    print len(E_TD), len(E_BO)
    

for e_td,e_bo,step in zip(E_TD,E_BO,steps):
    print>>of, step,e_td-e_bo
    
of.close()

plot="""set term pos color enhanced 'Helvetica' 20
        set output 'energy_above_BO.ps'
        set title 'Energy above BO surface'
        set key off
        HA=27.2114
        set yrange [0:]
        set xlabel 'output step'
        set ylabel 'E above BO (eV)'
        pl '%s' i 1 u 1:(($2)*HA) w l lw 4     
     """ %(ofile)
mix.gplot(plot,'energy_above_BO.gp')
viewer=os.environ.get('PSVIEWER')
mix.execute('%s energy_above_BO.ps' %viewer)


















#!/usr/bin/env python
"""
    Combine two wave function files from HOTBIT code into one wave function file.
    
    To be used in conjunction with TDTB calculations, where the initial
    wave function or state is taken as the eigenstate of two (or more) 
    separate systems.
    
    Usage:
        combine_wf.py input1 input2 output
        
        input1 - the first wave function file (wf.out from HOTBIT code)
        input2 - the second wave function file (-"-) 
        output - the name of the output file 
        
    Example:
        combine_wf.py wf1.out wf2.out wf.in
        
    The ordering of the wave function files is important. If you append the atoms
    of system 2 into the atoms.dat file after system 1, also the wave function
    files should be in the same order. The wave functions in the resulting file
    do not, on the other hand, have any definite order (order goes with the 
    occupation number, which can be quite random)
    The two systems should be separated so that they do not interact chemically,
    otherwise the orbitals will not be orthogonal.
    
    P. Koskinen 26.3 2007
"""
import box.mix as mix
import sys

try:
    in1,in2,out=sys.argv[1:]
except:
    print __doc__
    sys.exit(1)
    
n1  =int(mix.find_value(in1,'norb'))
n2  =int(mix.find_value(in2,'norb'))
occ1=mix.find_value(in1,'occ',fmt='matrix')
occ2=mix.find_value(in2,'occ',fmt='matrix')
ev1 =occ1[:,1]
occ1=occ1[:,0]
ev2 =occ2[:,1]
occ2=occ2[:,0]

of=open(out,'w')
norb=n1+n2
of.write('norb=%i\n' %norb)

eps=1E-5
i1=0
i2=0
occ=[]
for k in range(norb):
    if i1==n1: o1=0.0
    else:      o1=occ1[i1]
    if i2==n2: o2=0.0
    else:      o2=occ2[i2]
    if o1<eps and o2<eps:
        # do not write any wave functions if no occupations
        occ.append([0.0,0.0])
    elif o1>o2: 
        # state from system 1 is more occupied
        occ.append([o1,ev1[i1]])
        wf=mix.find_value( in1,'wf_'+str(i1+1),'strings' )
        of.write('wf_'+str(k+1)+'=\n')
        for i in range(len(wf)): 
            of.write( '%s\n' %(wf[i]) )
        of.write('\n')
        i1+=1
    else:
        # state from system 2 is more occupied. Add norb(system1) into basis state index.
        occ.append([o2,ev2[i2]])
        wf=mix.find_value( in2,'wf_'+str(i2+1),'strings' )
        of.write('wf_'+str(k+1)+'=\n')
        for i in range(len(wf)): 
            ind,w=wf[i].split()
            ind = int(ind)+n1
            of.write( '%i %s\n' %(ind,w) )
        of.write('\n')
        i2+=1  
    
of.write('occ=\n')
for k in range(norb): of.write('%f %f\n' %(occ[k][0],occ[k][1]) )
of.close()
    
    




    

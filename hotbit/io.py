#!/usr/bin/env python
"""
    Module for io-actions with the HOTBIT code.
"""
import box.mix as mix
import numpy as nu
vector=nu.array
err=mix.error_exit

def write_md_dat(mode=None,**kwargs):
    """ Make md.dat -file """
    
    defaults={'propagation':'BO',\
             'dt':1.0,\
             'mdsteps':1000,\
             'temp(K)':0.0,\
             'relax_time':0.0,\
             'e_temp':0.0005,\
             'quenching':'no',\
             'f_crit':5E-4,\
             'SCC':'no',\
             'SCC_crit':1E-3,\
             'SCC_mix':0.3}
             
    if mode=='black_box':
        defaults['mdsteps']=1000000000000000
        defaults['commmunicate']='yes'
        
    for key in kwargs:
        defaults[key]=kwargs[key]
           
    lines=""
    for key in defaults:
        lines=lines+key+'='+str(defaults[key])+'\n'
        
    f=open('md.dat','w')
    f.write(lines)
    f.close()
    
    
class BlackBox:
    """ Class for returning forces and energies for molecules. """
    
    def __init__(self,molecule=None,command='./source/hotbit'):
        """ Open the input and output pipes for the black box. """
        from os import popen2
        if molecule!=None:
            write_atoms_dat(molecule,'atoms.dat')
            md_dat({'communicate':'yes','mdsteps':10000})
        print 'Opening BlackBox pipe for communication...'
        self.inp,self.txt=popen2(command,'t')    
    
    def __del__(self):
        """ Close the pipe for the black box. """
        self.inp.write('1\n')
        self.inp.flush()
        print 'BlackBox pipe closed.'
    
    def get_energy_forces(self,molecule):
        """ Calculate forces and energy for given molecule. """
        N=molecule.get_N()
        lines=[]
        lines.append(str(N))
        r=molecule.get_positions()
        for pos in r:
            lines.append( mix.a2s(pos) )
        write_N_lines( self.inp,lines )
        ef = read_N_lines( self.txt )
        molecule.set_energy( float(ef[1]) )
        f=[]
        for line in ef[2:]:
            sl=line.split()
            f.append( vector([float(sl[0]),float(sl[1]),float(sl[2])]) )
        molecule.set_forces( f )
    
    def update(self,molecule):
        """ Update the molecule electronic structure. """
        self.get_energy_forces(molecule)
    
def read_N_lines(f):
    """ Read N lines from f, where N comes from the first line"""
    lines=[]
    lines.append( f.readline() )
    try:    
        N = int( lines[0] )
    except: 
        mix.error_exit("""line %s should be integer: number of lines in package.
                          Reading lines from file %s failed.
                       """ %(lines[0], f.name) )
    if N==1:        
        err('%s returned 1. HOTBIT returned error message.' %f.name)
    for i in range(N-1): 
        lines.append( f.readline() )
    return lines
    
    
def write_N_lines(f,lines):
    """ Write N lines into f, where N is written in the first line."""
    if len(lines)==0:
        mix.error_exit('[write_N_lines] Zero lines. Exit now.')
    N = int( lines[0] )
    for line in lines: 
        if line[-1:-2]!='\n': line = line+'\n'
        f.write(line) 
    f.flush()
    
def write_atoms_dat(molecule,file):
    """ Write the molecule into atoms.dat-file. """
    f=open(file,'w')
    f.write("nat=%i\n" %molecule.get_N())
    q=molecule.get_charge()
    f.write('extra_electrons=%f\n' %(-q))
    f.write("atoms=\n")
    for atom in molecule:
        print>>f, atom.get_symbol(), mix.a2s( atom.get_position() )
            
    minimum_set={'flag':'flags','v':'velocities','f':'forces'}
    for property in molecule.get_atom(0).list_properties():
        if property in minimum_set:
            name=minimum_set[property]
            print>>f, "\n\n",name,"="
            for atom in molecule:
                p=atom.get_property(property)
                if type(p)==type( vector([]) ): p=mix.a2s(p)
                print>>f, p
    f.close()
    
#def input_atoms_dat(self,file):
    #"""
    #Read molecule from atoms.dat-file.
    
    #* file -- the given atoms.dat -file
    #"""
    #fnd=mix.find_value
    #f=open(file,'r')
    #N=int(fnd(f,'nat'))
    #self.extra_electrons=float(fnd(f,'extra_electrons',default=0.0))
    #atoms=fnd(f,'atoms',fmt='strings')
    #for line in atoms:
        #el,x,y,z = line.split()
        #r=array([float(x),float(y),float(z)])
        #elem=Element(el,self.efile)
        #atom=Atom(elem,r)
        #self.add( atom )
        
    #if N!=self.N: error_exit('Error in atoms.dat file')
        
    #flags=fnd(f,'flags',fmt='strings',default='def')
    #if flags!='def': 
        #for i in range(self.N): self.atoms[i].act=int(flags[i])
    #else:            
        #for i in range(self.N): self.atoms[i].act=1
    #vels=fnd(f,'velocities',fmt='strings',default='def')
    #if vels!='def': 
        #for i in range(self.N):
            #self.atoms[i].v=array([float(x) for x in vels[i].split()]) 

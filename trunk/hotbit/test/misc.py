import os
from ase.data.molecules import molecule as ase_molecule
from ase import *

ddir='%s/hotbit/test/fixed_param/' %os.environ.get('HOTBIT_DIR')

default_elements={'C':ddir+'C.elm','Au':ddir+'Au.elm','H':ddir+'H.elm','N':ddir+'N.elm','Na':ddir+'Na.elm','O':ddir+'O.elm'}
default_tables={'AuAu':ddir+'Au_Au.par','CC':ddir+'C_C.par','HH':ddir+'H_H.par','NN':ddir+'N_N.par','CH':ddir+'C_H.par','CN':ddir+'C_N.par','NH':ddir+'N_H.par','NaNa':ddir+'Na_Na.par','OH':ddir+'O_H.par','OO':ddir+'O_O.par','CO':ddir+'C_O.par','AuC':ddir+'Au_C.par'}
default_param={'elements':default_elements,'tables':default_tables,'mixer':{'name':'Anderson','memory':3,'mixing_constant':0.2,'convergence':1E-7},'maxiter':1000,'width':0.02}


def molecule(name):
    try:
        return read('%s/hotbit/test/systems/%s.xyz' %(os.environ.get('HOTBIT_DIR'),name) )
    except:
        return ase_molecule(name)

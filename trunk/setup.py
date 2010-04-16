import os
import sys
import numpy as np

from distutils.core import setup, Extension
from distutils.sysconfig import get_config_var

sys.path += [ "." ]

from config import get_system_config

data_files = []
# data files & folders (folders nest only two levels down, not three anymore)
dirs = ['param','examples']
for dir in dirs:
   for item in os.listdir(dir):
       fullitem = os.path.join(dir,item)
       if os.path.isfile(fullitem):
           data_files.append((dir,[fullitem]))
       elif os.path.isdir(fullitem) and '.svn' not in fullitem:
           for item2 in os.listdir(fullitem):
               fullitem2 = os.path.join(fullitem,item2)
               if os.path.isfile(fullitem2):
                   data_files.append((fullitem,[fullitem2])) 


data_files.append(('.',['hotbit/hotbit']))

###
inc_dirs       = [ ]
lib_dirs       = [ ]
libs           = [ ]
extra_link     = [ ]
extra_compile  = [ ]

msgs = [ ]
get_system_config(inc_dirs,  libs, lib_dirs,
                  extra_link, extra_compile,
                  msgs)

# this is probably silly way of doing this:
version = '0.1'
revision=os.popen('svnversion .').readline()[:-1]
f=open('./hotbit/version.py','w').write('hotbit_version = "%s (svn=%s)"\n' %(version,revision))


setup(
    name         = "hotbit",
    url          = "https://trac.cc.jyu.fi/projects/hotbit",
    description  = "Density-functional tight-binding calculator for ASE",
    author_email = "pekka.koskinen@iki.fi",
    version      = version,
    packages     = [
        "box",
        "box.wedgeatoms",
        "hotbit",
        "hotbit.analysis",
        "hotbit.containers",
        "hotbit.parametrization",
        "hotbit.test"
        ],
    ext_modules  = [
        Extension(
            "_hotbit",
            [ "hotbit/c/_hotbit.c",
              "hotbit/c/geig.c", 
              "hotbit/c/slako.c",
              "hotbit/c/spherical.c",
              "hotbit/c/multipole.c" ],
            include_dirs        = inc_dirs,
            libraries           = libs,
            library_dirs        = lib_dirs,
            extra_compile_args  = extra_compile,
            extra_link_args     = extra_link
            )
        ],
    data_files = data_files    
    )

for msg in msgs:
    print msg

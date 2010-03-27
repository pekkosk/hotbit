import os
import sys

import numpy as np

from distutils.core import setup, Extension
from distutils.sysconfig import get_config_var

sys.path += [ "." ]

from config import get_system_config

# data files & folders
folders = ['param','param/inofficial','param/fixed_parameters','hotbit/test/systems',]
data_files = []
for folder in folders:
    files = []
    for file in os.listdir(folder):
        fullfile = os.path.join(folder,file)
        if os.path.isfile(fullfile):
            files.append(fullfile)
    data_files.append( (folder,files) )


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

setup(
    name         = "hotbit",
    url          = "https://trac.cc.jyu.fi/projects/hotbit",
    description  = "Density-functional tight-binding calculator for ASE",
    author_email = "pekka.koskinen@iki.fi",
    version      = "0.1",
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

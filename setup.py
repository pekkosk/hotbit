import os
import sys

import numpy as np

from distutils.core import setup, Extension
from distutils.sysconfig import get_config_var

sys.path += [ "." ]

from config import get_system_config

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
    version      = "0.1",
    packages     = [
        "box",
        "box.wedgeatoms",
        "hotbit",
        "hotbit.analysis",
        "hotbit.containers",
        "hotbit.parametrization",
        "hotbit.test",
        "hotbit.test.systems",
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
        ]
    )

for msg in msgs:
    print msg

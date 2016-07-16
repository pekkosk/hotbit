# Copyright (C) 2006 CSC-Scientific Computing Ltd.

# Please see the accompanying LICENSE file for further information.

# Adopted from GPAW

import os
import sys
from glob import glob
from os.path import join
from socket import gethostname 

# Always import numpy
import numpy


def find_file(arg, dir, files):
    #looks if the first element of the list arg is contained in the list files
    # and if so, appends dir to to arg. To be used with the os.path.walk
    if arg[0] in files:
        arg.append(dir)


def get_system_config(include_dirs, libraries, library_dirs,
                      extra_link_args, extra_compile_args, 
                      msg):

    include_dirs += [numpy.get_include()]
    include_dirs += ['c/libxc']

    machine = os.uname()[4]
    if machine == 'x86_64':

        #    _
        # \/|_||_    |_ |_|
        # /\|_||_| _ |_|  |
        #

        extra_compile_args += ['-Wall', '-std=c99']


        if 'MKL_ROOT' in os.environ:
            mklbasedir = [os.environ['MKL_ROOT']]
            libs = ['libmkl_intel_lp64.a']
            if mklbasedir != []:
                os.path.walk(mklbasedir[0],find_file, libs)
                libs.pop(0)
                if libs != []:
                    libs.sort()
                    libraries += ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_lapack', 'mkl_core',
                                  'iomp5', 'guide', 'mkl_def']#, 'mkl_def']
                    library_dirs += libs
                    msg +=  ['* Using MKL library: %s' % library_dirs[-1]]
        elif gethostname() in ['vuori1.csc.fi','vuori2.csc.fi']:
            libraries += ['acml']
            msg +=  ['* Using acml.']                      
        else:
            # Look for ACML libraries:
            acml = glob('/opt/acml*/g*64/lib')
            if len(acml) > 0:
                library_dirs += [acml[-1]]
                libraries += ['acml']
                if acml[-1].find('gfortran') != -1: libraries.append('gfortran')
                if acml[-1].find('gnu') != -1: libraries.append('g2c')
                extra_link_args += ['-Wl,-rpath=' + acml[-1]]
                msg += ['* Using ACML library']
            else:
                atlas = False
                for dir in ['/usr/lib', '/usr/local/lib', '/opt/lib', '/opt/local/lib']:
                    if glob(join(dir, 'libatlas.a')) != []:
                        atlas = True
                        break
                    if glob(join(dir, 'libatlas.dylib')) !=[]:
                        atlas = True
                        break
                if atlas:
                    libraries += ['lapack', 'atlas', 'blas']                    
                    library_dirs += [dir]
                    msg +=  ['* Using ATLAS library']
                else:
                    libraries += ['blas', 'lapack']
                    msg +=  ['* Using standard lapack']

    elif machine =='ia64':

        #  _  _
        # |_ |  o
        #  _||_||
        #

        extra_compile_args += ['-Wall', '-std=c99']
        libraries += ['mkl','mkl_lapack64']

    elif machine == 'i686':

        #      _
        # o|_ |_||_
        # ||_||_||_|
        #

        extra_compile_args += ['-Wall', '-std=c99']

        if 'MKL_ROOT' in os.environ:
            mklbasedir = [os.environ['MKL_ROOT']]
        else:
            mklbasedir = glob('/opt/intel/mkl*')
        libs = ['libmkl_ia32.a']
        if mklbasedir != []:
            os.path.walk(mklbasedir[0],find_file, libs)
        libs.pop(0)
        if libs != []:
            libs.sort()
            libraries += ['mkl_lapack',
                          'mkl_ia32', 'guide', 'pthread', 'mkl']#, 'mkl_def']
            library_dirs += libs
            msg +=  ['* Using MKL library: %s' % library_dirs[-1]]
            #extra_link_args += ['-Wl,-rpath=' + library_dirs[-1]]
        else:
            atlas = False
            for dir in ['/usr/lib', '/usr/local/lib']:
                if glob(join(dir, 'libatlas.a')) != []:
                    atlas = True
                    break
            if atlas:
                libraries += ['lapack', 'atlas', 'blas']
                library_dirs += [dir]
                msg +=  ['* Using ATLAS library']
            else:
                libraries += ['blas', 'lapack']
                msg +=  ['* Using standard lapack']

            # add libg2c if available
            g2c=False
            for dir in ['/usr/lib', '/usr/local/lib']:
                if glob(join(dir, 'libg2c.so')) != []:
                    g2c=True
                    break
                if glob(join(dir, 'libg2c.a')) != []:
                    g2c=True
                    break
            if g2c: libraries += ['g2c']

    elif sys.platform == 'darwin':

        extra_compile_args += ['-Wall', '-std=c99']
        include_dirs += ['/usr/include/malloc']

        extra_link_args += ['-framework Accelerate']

        if glob('/System/Library/Frameworks/vecLib.framework') != []:
            extra_link_args += ['-framework vecLib']
            msg += ['* Using vecLib']
        else:
            libraries += ['blas', 'lapack']
            msg +=  ['* Using standard lapack']

    return msg

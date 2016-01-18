import os
mklroot = os.environ['MKLROOT']

# include paths for compiling
#   inc_dirs = []

# library paths for linking
#   lib_dirs = []
lib_dirs = [mklroot + '/lib/intel64/', mklroot + '/../compiler/lib/intel64']

# libraries to link
#   libs = []
libs = ['mkl_intel_lp64',
        'mkl_intel_thread', 
        'mkl_lapack95_lp64', 
        'mkl_core',
        'mkl_def', 
        'iomp5']

# compiler options
#   extra_compile = []
extra_compile = ['-std=c99', '-O3']

# linker options
#   extra_link = []


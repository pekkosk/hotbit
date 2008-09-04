all:
	f2py -m eigensolver -c --verbose --f90exec='f95' --f90flags='-O3 -fPIC' -llapack -lblas eigensolver.f90 
	f2py -m slako -c --verbose --f90exec='f95' --f90flags='-O3 -fPIC' -llapack -lblas slako.f90
	f2py -m misc -c --verbose --f90exec='f95' --f90flags='-O3 -fPIC' -llapack -lblas slako.f90 misc.f90  


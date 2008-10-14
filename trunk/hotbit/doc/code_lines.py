import os
import pylab as pl
from box import mix
from time import time
folder=os.environ.get('HOTBIT_DIR')

python=0
fortran=0
for data in os.walk(folder):
    dirpath, dirnames, filenames = data
    if '.' in dirpath: continue
    for file in filenames:
        if file[-3:]=='.py':
            python+=len( open(dirpath+'/'+file).readlines() )
        if file[-2:]=='90':
            fortran+=len( open(dirpath+'/'+file).readlines() )
            
  

# append current lines into file
datafile='%s/hotbit/doc/code_lines.txt' %folder
o=open(datafile,'a')
t=time()
print>>o, t, python, fortran
o.close()

# plot the line info
t0=1223991504.5
data=mix.read(datafile)
data[:,0]-=t0
data[:,0]/=356*24*3600
xs, ys = pl.poly_between(data[:,0],0,data[:,1])
pl.fill(xs,ys,label='python')
xs, ys = pl.poly_between(data[:,0],0,data[:,2])
pl.fill(xs,ys,label='fortran')
pl.title('hotbit lines of code')
pl.xlabel('Years since 14.10 2008')
pl.ylabel('Lines of code')
pl.legend()
pl.savefig('code_lines.png')
pl.show()

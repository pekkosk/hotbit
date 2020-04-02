import os
import pylab as pl
import numpy as np   ## MS: needed for work-around of pylab.poly_between
from box import mix
from time import time,localtime
import datetime
folder='/home/pekkosk/workspace/hotbit' 
datafile=folder+'/hotbit/doc/code_lines.txt' 

python=0
fortran=0
c = 0
for data in os.walk(folder):
    dirpath, dirnames, filenames = data
    if '.' in dirpath: continue
    for file in filenames:
        if file[-3:]=='.py':
            python+=len( open(dirpath+'/'+file).readlines() )
        if file[-2:]=='90':
            fortran+=len( open(dirpath+'/'+file).readlines() )
        if file[-1:]=='c' or file[-1:]=='h':
            c+=len( open(dirpath+'/'+file).readlines() )
        

o=open(datafile,'a')
t=time()
print(t, python, fortran, c, file=o)
o.close()

fig = pl.figure()
ax = fig.add_subplot(111)

# plot the line info
t0=1223991504.5
data=mix.read(datafile)
timetuples = [localtime(sec) for sec in data[:,0]]
dates = [datetime.date(a[0],a[1],a[2]) for a in timetuples]
#data[:,0]-=t0
#data[:,0]/=356*24*3600
#pl.plot(dates,data[:,1])
#xs, ys = pl.poly_between(data[:,0],0,data[:,1])
## MS: incompatibility issue with matplotlib>=3.1
#xs, ys = pl.poly_between(dates,0,data[:,1])
#pl.fill(xs,ys,label='python')
#xs, ys = pl.poly_between(dates,data[:,1],data[:,1]+data[:,2])
#pl.fill(xs,ys,label='fortran')
#xs, ys = pl.poly_between(dates,data[:,1]+data[:,2],data[:,1]+data[:,2]+data[:,3])
#pl.fill(xs,ys,label='C')
pl.fill(np.append(dates,0), np.append(data[:,1],0), label='python')
pl.fill(np.append(dates,data[:,1]), np.append(data[:,1]+data[:,2],0), 
        label='fortran')
pl.fill(np.append(dates,data[:,1]+data[:,2]), 
        np.append(data[:,1]+data[:,2]+data[:,3],data[:,1]+data[:,2]),
        label='C')
pl.title('lines of code in hotbit')
pl.xlabel('Years since 14.10 2008')
pl.ylabel('Lines of code')
pl.legend(loc=2)

from matplotlib.dates import MonthLocator, WeekdayLocator, DateFormatter
from matplotlib.dates import MONDAY
mondays  = WeekdayLocator(MONDAY)
months   = MonthLocator() 
monthsFmt = DateFormatter('%b %y')

ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(monthsFmt)
ax.xaxis.set_minor_locator(mondays)

labels = ax.get_xticklabels()
pl.setp(labels,'rotation',45)

pl.savefig('code_lines.png')


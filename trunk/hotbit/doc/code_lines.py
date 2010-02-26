import os
import pylab as pl
from box import mix
from time import time,localtime
import datetime
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
xs, ys = pl.poly_between(dates,0,data[:,1])
pl.fill(xs,ys,label='python')
xs, ys = pl.poly_between(dates,0,data[:,2])
pl.fill(xs,ys,label='fortran')
pl.title('lines of code in hotbit')
pl.xlabel('Years since 14.10 2008')
pl.ylabel('Lines of code')
pl.legend()

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


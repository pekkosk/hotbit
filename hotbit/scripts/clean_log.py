#!/usr/bin/env python
"""
    Clean the hotbit.out -file from all the data from inside the main loop.
    (Leave only the initialization and final configuration data.)
        
    Usage:
        clean_log.py output
    
    P. Koskinen 12.4 2007
"""
import sys

try:
    out=sys.argv[1]
except:
    print __doc__
    sys.exit(1)
    
fi=open('hotbit.log','r')
lines=fi.readlines()
fi.close()

fo=open(out,'w')
flag=0
for line in lines:
    if line.find('Exit main iteration loop')>=0: flag=0
    if flag==0: fo.write(line)
    if line.find('Start main iteration loop')>=0: flag=1
    
fo.close()


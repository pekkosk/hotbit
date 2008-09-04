#!/usr/bin/env python
import os
from optparse import OptionParser
import box.mix as mix
clean=mix.clean_indentation
oeg=os.environ.get

email='pekka.koskinen@phys.jyu.fi'

# parse options
parser = OptionParser(usage='%prog [options] [script nnodes]')
#parser.add_option("-n", "--np", dest='np',help='number of processors',type='int',default=8)
parser.add_option("-m", "--memory", dest='memory',help='required memory in MB',type='int',default=0)
parser.add_option("-t", "--time", dest='time',help='wall time in seconds',default=344000)
opt, args = parser.parse_args()

# get environment
machine=os.environ.get('MACHINE') 
user=os.environ.get('USER')
if machine=='laptop' or machine=='nano110': machine='default'
if machine=='':
    mix.error_exit('Host not implemented.')

hotbit_exe=oeg('HOTBIT_EXE')
pythonpath=oeg('PYTHONPATH')
param=oeg('HOTBIT_PARAMETERS')


#if len(args)>1:
    #np=int(args[1])
#else:
    #np=opt.np                               #number of processors
np=1
cwd=os.path.abspath('.').split('/')[-1]     #current work dir= job name
name=args[0]                                #job name (the python script name)

of=open('run_%s' %name,'w')

if machine=='murska':
    t=mix.sec_to_time(opt.time)  #time limit ("dd:hh:mm:ss")
    time='%02i:%02i' %(t[0]*24+t[1],t[2])
    #if opt.memory==0:
        #mix.error_exit('Required memory must be specified')
    #mem_core=int(opt.memory/np)*1000 
    batch=clean("""
        #!/bin/bash
        #BSUB -n %(np)i
        #BSUB -W %(time)s
        #BSUB -J %(cwd)s
        #BSUB -u %(email)s
        #BSUB -N
        #BSUB -e %(name)s_%%J.err
        #BSUB -o %(name)s_%%J.out
        
        module load gpaw
        export HOTBIT_EXE="%(gpaw_setup_path)s"
        export PYTHONPATH="%(pythonpath)s"
        export HOTBIT_PARAMETERS="%(param)s"
        python %(name)s """ %vars())
        
    

elif machine=='batman':
        time=opt.time
        batch=clean("""
            #PBS -N %(name)s
            #PBS -l ncpus=%(np)i
            #PBS -l walltime=%(time)s
            #PBS -m bea
            #PBS -M %(email)s
            #PBS -e %(name)s.err ##PBS_JOBID?
            #PBS -o %(name)s.out
            
            # change to the directory where you submitted the job
            cd %(cwd)s
                        
            export HOTBIT_EXE="%(hotbit_exe)s"
            export PYTHONPATH="%(pythonpath)s"
            export HOTBIT_PARAMETERS="%(param)s"
            
            python %(name)s >> %(name)s_%%J.out
            exit 0""" %vars())
print>>of, batch
of.close()

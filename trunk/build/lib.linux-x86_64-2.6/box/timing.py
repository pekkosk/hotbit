# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

from time import time, asctime
import numpy as np

class Timer:
    """ Class for making timing in nested fashion. 
    
    txt=open(...)
    tm=Timer('main program',txt)
    ...
    tm.start('sub')
    tm.stop('sub')
    ...
    tm.start('sub2')
    tm.start('sub2-sub')
    tm.stop('sub2-sub')
    tm.stop('sub2')
    ...
    tm.summary()
    
    Timing can be nested, but still strictly hierarchial. 
    (E.g. if function (with timer) is called from different places (with different timers),
    timers will overlap and an error will result.)    
    Cannot be used to recursive functions.
    """
    
    def __init__(self,label,level=0,txt=None,enabled=True):
        """ Init process timer with given label and output. 
        
        Parameters:
        -----------
        label:      label for process (e.g. 'integration')
        
         
       
        
        level refers to deepness of nesting. 
        
        If not enabled, disable all timing stuff (start,stop,summary) for speeding up.        
        """
        self.label=label
        self.timers=[]
        self.first=time()
        if txt==None:
            from sys import stdout
            self.txt=stdout
        else:
            self.txt=txt            
        self.level=level
        self.running=False
        self.enabled=enabled
        if level==0:
            self.start()
        #self.durations=[]
        self.elapsed_time=0.0
        self.nr_of_calls=0
        self.smry=False

            
    def is_running(self):
        return self.running 
           
    def get_time(self):
        if self.running:
            self.stop()          
        return self.elapsed_time
        
    def get_calls(self):
        if self.running:
            self.stop()
        return self.nr_of_calls        
        
    def get_timer(self,label):
        """ Does timer or any of it's timers have given label (recursively)? """
        if self.label==label:
            return self               # this has given label
        elif len(self.timers)==0:
            return None                # no given label, no child timers
        else:
            for timer in self.timers:
                tr=timer.get_timer(label)
                if tr!=None: 
                    return tr
            return None 
                  
    def get_level(self):
        return self.level            
                  
    def get_outmost_running(self):
        """ Return the most nested running timer. 
        If timer is running and has no child timers, itself is the outmost timer.
        If if no child timers are running, return itself.
        """
        if not self.running:
            raise AssertionError('Timer not running while searching for outmost running timer.')
        if len(self.timers)==0:
            return self
        else:
            for timer in self.timers:
                if timer.is_running():
                    return timer.get_outmost_running()
            return self        
                    
    def add_subtimer(self,tr):
        self.timers.append(tr)                    
                        
    def start(self,label=None):
        """ Start timer itself or child timer with given label. """   
        if not self.enabled: return
        if label==None:
            if self.running:
                raise AssertionError('Timer %s already running!' %self.label)
            self.t1=time()
            self.running=True
        else:            
            tr=self.get_timer(label)   
            if not self.running:
                 raise AssertionError('Timer %s cannot make running child timers; itself is not running!' %self.label)
            if tr==None:
                outmost=self.get_outmost_running()
                tr=Timer(label,txt=self.txt,level=outmost.get_level()+1)    
                outmost.add_subtimer(tr)
            if tr.is_running():
                raise AssertionError('Timer %s is already running!' %label)                
            tr.start()
            
    def stop(self,label=None):
        if not self.enabled: return
        if label==None:
            if not self.running:
                raise AssertionError('Timer %s cannot be stopped; it is not running!' %self.label)
            self.running=False
            self.t2=time()
            self.elapsed_time += self.t2-self.t1
            self.nr_of_calls += 1
        else:
            tr=self.get_timer(label)
            if tr==None:
                raise AssertionError('Timer %s cannot be stopped; it does not exist!' %label)
            tr.stop()            
        
    def get_summary(self,total,partial,dict):    
        dt, calls=self.get_time(), self.get_calls()
        dt_sub=0.0
        txt_sub=''
        subs=False
        for timer in self.timers:
            dt2, txt2, dict=timer.get_summary(total,dt,dict)
            txt_sub+=txt2
            dt_sub+=dt2
            subs=True
        sub_covered=dt_sub/dt*100
        
        procent1=dt/total*100.0
        x=int(np.round(procent1*0.3))
        bar='|'+str(self.level)*x+'|'
        procent2=dt/partial*100.0
        txt='../'*self.level+'%-20s' %(self.label) 
        if subs:
            txt+=' '*(10-3*self.level)+'%12.3f %9i (%5.1f %%,%5.1f %%) %5.1f %% %s\n' %(dt,calls,procent2,sub_covered,procent1,bar)
        else:
            txt+=' '*(10-3*self.level)+'%12.3f %9i (%5.1f %%        ) %5.1f %% %s\n' %(dt,calls,procent2,procent1,bar)
        dict[self.label]=dt
        
        txt+=txt_sub            
        return dt, txt, dict
        
    def get_timings(self):
        """ Return dictionary of processes and their execution times. """       
        if not self.enabled: return
        if not self.smry:
            raise AssertionError('Summary must be performed before get_timing.')
        return self.dict
        
        
    def summary(self):  
        if not self.enabled: return  
        self.last=time()
        total=self.last-self.first
        dict={}
        dt, txt, self.dict=self.get_summary(total,total,dict)
                
        print>>self.txt, '\nTiming:'
        print>>self.txt, '            label                    time     calls    %sub  %covered   %tot'
        print>>self.txt, '-'*79
        print>>self.txt, txt,
        print>>self.txt, '-'*79
        print>>self.txt, 'total time %12.3f seconds      %s' % (total, self.human_readable_time(total))
        print>>self.txt, asctime()
        self.smry=True
        self.txt.flush()


    def human_readable_time(self, seconds):
        seconds = int(round(seconds))
        hours = seconds/3600
        seconds -= hours*3600
        mins = seconds/60
        seconds -= mins*60
        return "%i h  %i min  %i sec" % (hours, mins, seconds)


class OneTimer:
    def __init__(self,key):
        self.key=key
        self.t0=time.time()   
        self.durations=[]
        self.running=False
        
    def start(self):
        self.t1=time.time()
        self.running=True
       
    def stop(self):
        assert self.running
        self.t2=time.time()
        self.durations.append( self.t2-self.t1 )
        self.running=False
        
    def is_running(self):
        return self.running
        
    def get_time(self):
        if self.running==True:
            self.stop()
        return sum(self.durations)
    
    def get_calls(self):
        return len(self.durations)
        

          
            

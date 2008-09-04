from time import time, asctime
import numpy as nu

class Timer:
    """ Class for making timing in nested fashion. 
    
    out=open(...)
    tm=Timer('main program',out)
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
    
    def __init__(self,label,out,level=0,enabled=True):
        """ Init process timer with given label and output. level refers to deepness of nesting. 
        
        If not enabled, disable all timing stuff (start,stop,summary) for speeding up.        
        """
        self.label=label
        self.timers=[]
        self.first=time()
        self.out=out
        self.level=level
        self.running=False
        self.enabled=enabled
        if level==0:
            self.start()
        self.durations=[]
        self.smry=False

            
    def is_running(self):
        return self.running 
           
    def get_time(self):
        if self.running:
            self.stop()          
        return sum(self.durations)
        
    def get_calls(self):
        if self.running:
            self.stop()
        return len(self.durations)        
        
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
                tr=Timer(label,self.out,outmost.get_level()+1)    
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
            self.durations.append( self.t2-self.t1 )
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
        x=int(nu.round(procent1*0.3))
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
                
        print>>self.out, '\nTiming:'
        print>>self.out, '            label                    time     calls    %sub  %covered   %tot'
        print>>self.out, '-'*79
        print>>self.out, txt,
        print>>self.out, '-'*79
        print>>self.out, 'total time %12.3f' %total
        print>>self.out, asctime()
        self.smry=True
        self.out.flush()
                


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
        

#class Timer:
    #def __init__(self,label,out):
        #self.label=label
        #self.procs={}
        #self.t0=time.time()
        #self.out=out
    
    #def start(self,key):
        #""" Start process with given key. """        
        #if key not in self.procs:
            #self.procs[key]=[OneTimer(key),{}]
        #for k2 in self.procs:
            #assert not self.procs[key].is_running()
        #self.procs[key].start()
        
    #def stop(self,key):
        #""" Stop process with given key. """
        #self.procs[key].stop()
        
    #def __del__(self):
        #self.summary()
        
    #def summary(self):
        #self.t1=time.time()
        #total=self.t1-self.t0        
            
        #print>>self.out, '\nTiming:'
        #print>>self.out, '-'*79
        #for key in self.procs:
            #dt=self.procs[key].get_time()
            #calls=self.procs[key].get_calls()
            #procent=dt/total*100.0
            #x=int(procent*0.3)
            #bar='|'+'-'*x+'|'
            #print>>self.out, '%25s %12.3f (%5.1f %%) %5i %s' %(key,dt,procent,calls,bar)
        #print>>self.out, '-'*79
        #print>>self.out, 'total time %12.3f' %total
        #print>>self.out, time.asctime()
            
            

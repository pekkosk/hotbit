#!/usr/bin/env python
import box.mix as mix
from numpy import *

class Gnuplot:
        
    def __init__(self,gfile,ofile,r,c,title='',xlabel='',ylabel=''):
        self.gfile=gfile
        self.ofile=ofile
        self.tfile=gfile+'_tmp'
        self.tmp=open(self.tfile,'w')
        self.c=c
        self.r=r
        self.N=r*c
        self.dy=1./r
        self.dx=1./c
        if self.N==1: 
            self.multi=False
        else:
            self.multi=True
            
        self.start="""
            set term pos col portrait enhanced
            set output "%(ofile)s"
            set title  "%(title)s"
            set xlabel "%(xlabel)s"
            set ylabel "%(ylabel)s"
            set tmargin 0
            set bmargin 3
            set rmargin 0
            set lmargin 4
            set size 1.,1.
        """ %vars()
        if self.multi: self.start+="set multiplot \n"
        self.plots=[ ['\n'] for x in range(r*c)]
        self.i=-1
        self.flags=[0 for x in range(self.N)]
    
    def plot(self,r,c,d,title='',key=''):
        mix.write(self.tmp,d,'%20.10g')
        self.tmp.write('\n\n')
        self.i+=1
        n=(c-1)*self.r+r-1
        
        if self.flags[n-1]==0:
            x=(c-1)*self.dx
            y=1.0-r*self.dy
            self.plots[n].append( "set origin %f,%f \n" %(x,y) )
            self.plots[n].append( "set size %f,%f \n" %(self.dx,self.dy) )
            self.plots[n].append( "set title '%s'\n" %title )       
            cmd='pl'
        else:
            cmd='repl'
            
        line="%s '%s' i %i u 1:2 w l t '%s'\n" %(cmd,self.tfile,self.i,key)
        self.plots[n].append(line)
        self.flags[n-1]=1 
        
    def execute(self,remove=True):
        self.tmp.close()
        
        # make the gnuplot file ...
        f=open(self.gfile,'w')
        f.writelines(self.start)
        for lines in self.plots:
            f.writelines( lines )
        f.close()
        
        # ... and call gnuplot
        mix.execute('gnuplot %s' %self.gfile)
        from os import remove
        #if remove: remove(self.tfile)
        
        
if __name__=='__main__':
    print 'test gnuplot'
    gp=Gnuplot('plot.gp','plot.ps',5,2,\
               title='R_{nl}(r)',xlabel='r (a_B)',ylabel='R_{nl}(r)')
    x=array([[1,2,3,4],[1,3,4,2]]).transpose()
    gp.plot(1,1,x,'koe')
    
    x=array([[1,2,3,4],[1,3,4,2]]).transpose()
    gp.plot(2,1,x,'koe')
    
    x=array([[1,2,3,4],[2,3,4,2]]).transpose()
    gp.plot(1,2,x,'koe')
    
    x=array([[1,2,3,4],[3,3,4,1]]).transpose()
    gp.plot(3,2,x,'koe')
    gp.execute()
            
            
            
            
            
     
            
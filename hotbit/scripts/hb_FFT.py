#!/usr/bin/env python
"""
    Calculate the Fourier spectrum from dipole oscillations.
    
    Usage:
        
        hb_FFT.py comp emax [k]
        
        * comp - the component of dipole moment which is used in transformation.
        
        * emax - max energy in eV to be plotted
        
        * k -- (optional) The decay time (in units eV/hbar) to reduce the errors
          from finite sampling. E.g. '0.3' eV/hbar
        
    Example:
        
        hb_FFT.py x 4.0
        
    Pekka Koskinen 21.3 2007
"""
import box.mix as mix
import os,sys
from numpy import *
from FFT import *
print 'Calculating the Fourier-transformed signal of dipole oscillations...'

try:
    component=sys.argv[1]
    emax=float(sys.argv[2])
except:
    print __doc__
    sys.exit(0)

try:
    k=float(sys.argv[3])
except:
    k=0.0
    
hbar=0.02342
k=k/27.2114/hbar
tau=1/k
if k>1E-10: print 'damping time (a.u.)',tau
M=20
print 'average over',M,'damped FFTs'

components={'x':1,'y':2,'z':3}
print "Make FFT using dipole component",component
m   = mix.read('dipolemom.out')
x   = m[:,0]
n   = size(x)
y   = m[:,components[component]]
x0_list  = linspace( x[0],x[-1],M )
if max(x)-min(x)<3*tau: mix.error_exit('1/k too large compared to the whole time-span')

print 'weight f(t) by exp(-k|t-t0|),t0=%12.4f ...%12.4f' %(x0_list[0],x0_list[-1])
sm  = 0.0*x
for x0 in x0_list:
    
    y   = y*exp( -k*abs(x-x0) )
    t   = abs(fft(y))
    sm += t
dt  = (x[1]-x[0])

# total time span T=n*dt
# omega_i = 2*pi/(T/i)=i*2*pi/(n*dt)
# energy (eV)    = hbar*omega_i*27.2114
# energy (cm^-1) = i/dt*32296
e = arange(0,n)*2*pi/(n*dt) * 0.0234218*27.21139


print 'output of the Fourier-transformed signal'
f=open('fft.out','w')
for i in range(1,n/2+2):
    f.write('%12.6f %12.6f\n' %(e[i],sm[i]) )
f.close()

file='plot_fft'
plot="""set term pos color enhanced "Helvetica" 20
    set output '%s.ps'
    set title 'FFT of the dipole signal'
    set xlabel 'energy (eV)'
    set ylabel 'amplitude (arb.units)'
    set xrange [0:%f]
    set yrange [0:]
    set key off
    pl 'fft.out' u 1:2 w l lw 4 """ %(file,emax)
mix.gplot(plot,file+'.gp')
psviewer=os.environ.get('PSVIEWER')
mix.execute('%s %s.ps &' %(psviewer,file) )


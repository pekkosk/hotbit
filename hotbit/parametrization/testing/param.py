from hotbit.parametrization import KSAllElectron
from hotbit.parametrization.util import IP_EA
from box.mix import Timer
#import cgitb; cgitb.enable()
import pylab as pl

#IP=ionization_potential('C',remove='2p')
IP, EA=IP_EA('O',remove_orb='2p',add_orb='2p',remove=1.0,add=0.5)
print(IP,IP*27.2114,IP*27.2114/0.01036)
print(EA,EA*27.2114,EA*27.2114/0.01036)
raise SystemExit

for i,element in enumerate(['H','C','Ti','Au']):
    e={}
    per=[100,150,200,250,300,400,500,1000]
    for pernode in per:
        atom=KSAllElectron(element,pernode=pernode)
        v=atom.get_valence()
        if e=={}:
            for val in v:
                e[val]=[]
        atom.run()
        for val in v:
            e[val].append(atom.get_eigenvalue(val))
        
    pl.subplot(2,2,i+1) 
    for val in v:
        #pl.plot(per,e[val]-e[val][-1],label=val)
        pl.semilogy(per,e[val]-e[val][-1]+1E-6,label=val)
    pl.legend()
    pl.xlabel('points/node')
  
#pl.savefig('grid_convergence.png')
pl.show()
    

#atom.plot_Rnl(screen=True)
#atom.write_functions('H.elm')
#atom.plot_Rnl(screen=True)
#print 'valence orb', atom.get_valence()
#print 'valence', atom.get_valence_energies()
#print 'symb',atom.get_symbol()
#print 'comment',atom.get_comment()
tm()





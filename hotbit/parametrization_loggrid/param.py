from hotbit.parametrization import KSAllElectron
from box.mix import Timer
#import cgitb; cgitb.enable()

tm=Timer()
#atom=KSAllElectron('H',etol=1E-6,occu={'1s':1},convergence=1E-4)
atom=KSAllElectron('H',verbose=True,etol=1E-10,convergence=1E-5)
#atom.solve_eigenstates()
atom.solve_ground_state()
atom.write_functions('H.elm')
#atom.plot_Rnl(screen=True)
#print 'valence orb', atom.get_valence()
#print 'valence', atom.get_valence_energies()
#print 'symb',atom.get_symbol()
#print 'comment',atom.get_comment()
tm()





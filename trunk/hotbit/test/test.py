import os
from time import time
import os

tests=['atom_dimer.py','standard_set.py','forces.py','external_field.py','parametrization.py',\
       'Au_chain.py','graphene.py','CNT5_0_chiral.py','polyethene_twisted.py','C6H6_wedge.py',\
       'linear_response.py','C6H6_cell_game.py','save_load.py','copy_calculator.py',\
       'mulliken.py','multipole_operations.py','periodicity.py']

       
skip = ['save_load.py','copy_calculator.py']
start = time()

pth=os.environ.get('HOTBIT_DIR')

for test in tests:
    if test in skip:
        print 'test', test,'skipped...'
        continue 
    try:
        file = os.path.join(pth,'hotbit','test',test)
        t1 = time()
        ret=os.system('python %s' %file)
        elapsed = time()-t1
        if ret!=0:
            print test,'returned',ret,'and FAILED!'
        else:
            print '%-25s OK. (%.1f seconds)' %(test,elapsed)
    except:
        print test,'ERROR!'

stop = time()
print "Total time elapsed: %.0f seconds." %(stop-start)

import os
import traceback
import subprocess
from time import time

tests = [
    'atom_dimer.py',
    'standard_set.py',
    'forces.py',
    'external_field.py',
    'parametrization.py',
    'Au_chain.py',
    'graphene.py',
    'CNT5_0_chiral.py',
    'polyethene_twisted.py',
    'C6H6_wedge.py',
    'linear_response.py',
    'C6H6_cell_game.py',
    'mulliken.py',
    'd-orbital_rotation_yaxis.py',
    'multipole_operations.py',
    'periodicity.py',
    'madelung_constants.py',
    'mio.py']

       
skip = []
add_env = {
    'mio.py': 'MIO_0_1'
    }
    
start = time()

pth=os.environ.get('HOTBIT_DIR')

for test in tests:
    if test in skip:
        print('test', test,'skipped...')
        continue 
    if test in add_env:
        if not add_env[test] in os.environ:
            print('test', test, 'requires environment variable', \
                add_env[test], 'skipped...')
            continue
    try:
        if pth is None:
            file = test
        else:
            file = os.path.join(pth,'hotbit','test',test)
        t1 = time()
        ret=os.system('python3 '+file)
        elapsed = time()-t1
        if ret!=0:
            print(test,'returned',ret,'and FAILED!')
        else:
            print('%-25s OK. (%.1f seconds)' %(test,elapsed))
    except:
        print(test,'ERROR!')
        traceback.print_exc()

stop = time()
print("Total time elapsed: %.0f seconds." %(stop-start))

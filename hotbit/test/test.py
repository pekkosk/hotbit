import os
tests=['atom_dimer.py','standard_set.py','forces.py','parametrization.py','linear_response.py','save_load.py', 'copy_calculator.py']

for test in tests:
    try:
        ret=os.system('python %s' %test)
        if ret!=0:
            print test,'returned',ret,'and FAILED!'
        else:
            print test,'OK.'
    except:
        print test,'ERROR!'

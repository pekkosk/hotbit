from hotbit.fortran.slako import fast_slako_transformations as fast_slako_transformations_for
from _hotbit import fast_slako_transformations as fast_slako_transformations_c

from math import sqrt

import numpy as nu

###

for i in range(100):
    rhat  = nu.random.random( 3 )
    rhat /= sqrt(nu.dot(rhat, rhat))

    h = nu.random.random( 14 )
    s = nu.random.random( 14 )
    dh = nu.random.random( ( 14, 3 ) )
    ds = nu.random.random( ( 14, 3 ) )

    ht1, st1, dht1, dst1 = fast_slako_transformations_for(rhat, 1.0, 9, 9, h, s, dh, ds)
    ht2, st2, dht2, dst2 = fast_slako_transformations_c(rhat, 1.0, 9, 9, h, s, dh, ds)

    eht = nu.max(nu.abs(ht1-ht2))
    est = nu.max(nu.abs(st1-st2))
    edht = nu.max(nu.abs(dht1-dht2))
    edst = nu.max(nu.abs(dst1-dst2))

    if eht > 1e-12 or est > 1e-12 or edht > 1e-12 or edst > 1e-12:
        print "error(ht) = %e" % eht
        print "error(st) = %e" % est
        print "error(dht) = %e" % edht
        print "error(dst) = %e" % edst


    


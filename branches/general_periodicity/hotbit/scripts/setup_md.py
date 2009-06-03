#!/usr/bin/env python
"""
    Make md.dat with given scheme for HOTBIT simulations.
    
    Usage:
        setup_md.py {scheme} {output}
        
        scheme - the simulation scheme. Available schemes are:
            MD                  - const T MD simulation
            quench              - optimize given structure with SCC
            linear_response     - calculate linear optical response
            response_fft        - calculate response for small laser pulse
            once                - one SCC calculation without propagation
            laser               - TDTB calculation with laser
            communicate         - use HOTBIT as a black box
     
    Example:
        setyp_md.py quench md.dat
        
    P. Koskinen 16.4 2007
"""        
import sys
import box.mix as mix

try:
    scheme,output=sys.argv[1:]
except:
    print __doc__
    sys.exit(1)
if scheme=='MD':
    out="""propagation=BO
    dt          = 1
    mdsteps     = 1000
    temp(K)     = 100
    relax_time  = 200
    e_temp      = 0.0007
    quenching   = no
    SCC         = yes
    SCC_crit    = 1E-4
    SCC_mix     = 0.3
    """
elif scheme=='communicate':
    out="""propagation=BO
    communicate = yes
    dt          = 0
    mdsteps     = 1
    temp(K)     = 0
    relax_time  = 0
    e_temp      = 0.0007
    quenching   = no
    SCC         = yes
    SCC_crit    = 1E-4
    SCC_mix     = 0.3
    """
elif scheme=='linear_response':
    out="""propagation=BO
    dt          = 0
    mdsteps     = 1
    temp(K)     = 0
    relax_time  = 0
    e_temp      = 0.0005
    quenching   = no
    SCC         = yes
    SCC_crit    = 1E-4
    SCC_mix     = 0.3
    out_optical = yes
    optical_de  = 0.5 # beware of this!
    calculate_forces=no
    #run hb_optical.py after calculation
    """
elif scheme=='quench':
    out="""propagation=BO
    dt          = 1
    mdsteps     = 1000
    temp(K)     = 0
    relax_time  = 0
    e_temp      = 0.0006
    quenching   = yes
    SCC         = yes
    SCC_crit    = 1E-4
    SCC_mix     = 0.3
    """
elif scheme=='once':
    out="""propagation=BO
    dt          = 0
    mdsteps     = 1
    temp(K)     = 0
    relax_time  = 0
    e_temp      = 0.0006
    quenching   = no
    SCC         = yes
    SCC_crit    = 1E-4
    SCC_mix     = 0.3
    out_charges=yes
    """
elif scheme=='laser':
    out="""propagation=TDTB
    dt          = 0.03
    mdsteps     = 1000
    temp(K)     = 0
    relax_time  = 0
    e_temp      = 0.0006
    quenching   = no
    SCC         = yes
    SCC_crit    = 1E-7
    SCC_mix     = 0.3
    
    laser=yes
    laser_hw(eV)=2.0
    laser_t0=2
    laser_dt=4
    laser_E =0 0 1
    laser_flux(W/cm2)=1E10
    
    laser2=no
    laser2_hw(eV)=2.0
    laser2_t0=10
    laser2_dt=4
    laser2_E =0 0 1
    laser2_flux(W/cm2)=1E10
    
    out_dipolemom=yes
    out_charges=yes
    """
elif scheme=='response_fft':
    out="""propagation=TDTB
    dt          = 0.03
    mdsteps     = 10000
    temp(K)     = 0
    relax_time  = 0
    e_temp      = 0.0006
    quenching   = no
    SCC         = yes
    SCC_crit    = 1E-7
    SCC_mix     = 0.3
    laser=yes
    laser_hw(eV)=2.0
    laser_t0=0.05
    laser_dt=0.1
    laser_E =1 1 1
    laser_flux(W/cm2)=1E11
    out_dipolemom=yes
    out_charges=yes
    calculate_forces=no
    """
else: mix.error_Exit('Scheme not defined.')
    
o=open(output,'w')
for line in out.splitlines():
    o.write(line.lstrip()+'\n')
o.close()
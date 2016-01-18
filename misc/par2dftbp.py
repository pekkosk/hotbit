#
#	Script to transform .par to .skf files to be used for DFTB+.
#
from helper.mix import find_value
from numpy import *
from helper.interpolation import SplineFunction
from box.mix import fit
from pylab import *
from scipy.optimize import fmin


def vrep_poly(r,r0,c):
    """ Repulsive potential as polynomial.

    Vrep(r) = sum_i=0^n c_i (r0-r)**n
    """
    n = len(c)
    assert all(r>=r0)
    return array([c[i]*(r-r0)**i for i in range(n)]).sum(axis=0)

def vrep_to_spline_coefficients(rep):
    """ 
        Transform simple grid data to spline form used by DFTB+

    Input:
    rep:    rep[:,0]=r and rep[:,1]=v(r)
    
    Output:
    spline: List with spline data
            Each line with: start end c0 c1 c2 c3
            where v is interpolated between [start,end] using
            v(r) = c0 + c1*(r-start)+ c2*(r-start)**2 + c3*(r-start)**3
            The last line has coefficients up to c5.
    """
    print rep.shape
    v = SplineFunction(rep[:,0],rep[:,1])
    n = rep.shape[0]
    m = 6
    spline = []
    for i in range(n-1):
        r0, r1 = rep[i,0], rep[i+1,0]
        v0, v1 = rep[i,1], rep[i+1,1]
        rlist = linspace(r0,r1,m)
        vlist = v(rlist)
        if i==n-2:
            nc = 1+5
        else:
            nc = 1+3
        pguess = zeros(nc)
        pguess[0] = 0.5*(v0+v1)
        pguess[1] = (v1-v0)/(r1-r0)
        def chi2(p,x,y):
            err = y-vrep_poly(x,r0,p)
            return sum(err**2)
        p = fmin(chi2,pguess,args=(rlist,vlist),xtol=1E-12,disp=False)
        spline.append([r0,r1]+list(p))
    return spline


def read_par(el1,el2,filename):
    """ Read parameter files. """
    f = open(filename,'r')
    tables = {}
    t12 = find_value(f,'%s_%s_table' %(el1,el2),fmt='matrix')
    t21 = find_value(f,'%s_%s_table' %(el2,el1),fmt='matrix')
    tables['N'] = len(t12[:,0])
    dr = t12[1,0]-t12[0,0]
    dr2 = t21[1,0]-t21[0,0]
    assert abs(dr-dr2)<1E-13 
    tables['dr'] = dr
    tables['grid'] = t12[:,0]
    tables['%s%s' %(el1,el2)] = t12[:,1:] 
    tables['%s%s' %(el2,el1)] = t21[:,1:]

    rep = find_value(f,'repulsion',fmt='matrix')
    return tables, rep


def read_elm(filename,orbitals):
    """ Read element data form .elm-file. """
    f = open(filename,'r')
    data = {}
    data['symbol'] = find_value(f,'symbol')
    data['energies'] = [float(find_value(f,'epsilon_%s' %orb)) for orb in orbitals]
    data['U'] = float(find_value(f,'U'))
    data['FWHM'] = float(find_value(f,'FWHM'))
    f.close()
    return data    


def write_skf(el1,el2,data1,data2,tables,rep,spline,dr=0.1):
    """ Write .skf files for el1-el2 interactions. 

    Input:
    el1:    element 1 symbol
    el2:    element 2 symbol 
    data1:  dictionary for element 1 data
    data2:  dictionary for element 2 data
    tables: dictionary for Slater-Koster table data
    rep:    repulsive potential on a grid
    spline: spline in appropriate format already
    dr:     output grid spacing
    """
    if el1==el2:
        filename = '%s-%s.skf' %(el1,el2)
        o = open(filename,'w')
        ng = int(ceil(tables['grid'][-1]/dr))
        print>>o, '%.15f, %i, ' %(dr,ng)
        e = data1['energies']
        u = data1['U']
        occ = data1['occupations']
        # Ed Ep Es SPE Ud Up Us fd fp fs (orbital energies and Hubbard U)
        print>>o, e[0], e[1], e[2], 0., u, u, u, occ[0], occ[1], occ[2]
        # mass c2 c3 c4 c5 c6 c7 c8 c9 rcut d1 d2 d3 d4 d5 d6 d7 d8 d9 d10
        # where c* are vrep coefficients (unless spline) and d* are not used
        print>>o, data1['mass'], 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
        #
        # interpolate SlaKo tables
        #
        tab = tables['%s%s' %(el1,el2)] #[i,j]
        f = [SplineFunction(tables['grid'],tab[:,i]) for i in range(20)]
        for i in range(ng):
            for j in range(20):
                print>>o, f[j]((i+1)*dr),
            print>>o        
        #
        # repulsion
        #
        print>>o, 'Spline'
        print>>o, len(spline), spline[-1][1] #nInt, cutoff
        print>>o, 0,-1E19,spline[0][2]
        for sp in spline:
            print>>o, ' '.join([str(x) for x in sp])           
        print>>o
        o.close()

    else:
        for e1,e2 in [(el1,el2),(el2,el1)]:
            filename = '%s-%s.skf' %(e1,e2)
            o = open(filename,'w')
            ng = int(ceil(tables['grid'][-1]/dr))
            print>>o, '%.15f, %i, ' %(dr,ng)
            print>>o, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.

            #
            # interpolate SlaKo tables
            #
            tab = tables['%s%s' %(e1,e2)] #[i,j]
            f = [SplineFunction(tables['grid'],tab[:,i]) for i in range(20)]
            for i in range(ng):
                for j in range(20):
                    print>>o, f[j]((i+1)*dr),
                print>>o        
            #
            # repulsion
            #
            print>>o, 'Spline'
            print>>o, len(spline), spline[-1][1] #nInt, cutoff
            print>>o, 0,-1E19,spline[0][2]
            for sp in spline:
                print>>o, ' '.join([str(x) for x in sp])           
            print>>o
            o.close()




el1, elm1, orb1 = 'Au', 'Au.elm', ['5d','6p','6s']
el2, elm2, orb2 = 'C', 'C.elm', ['2p','2s']
parfile = 'Au_C.par'

data1 = read_elm(elm1,orb1)
data2 = read_elm(elm2,orb2)
tables, rep = read_par(el1,el2,parfile)

data1['occupations'] = [10.,0.0,1.0]
data1['mass'] = 196.97
data2['occupations'] = [2.0,2.0]
data2['mass'] = 12.01

spline = vrep_to_spline_coefficients(rep)

write_skf(el1,el2,data1,data2,tables,rep,spline)


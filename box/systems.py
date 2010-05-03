#"""
    #Make various atomic structures into xyz-files.
    
    #P. Koskinen 31.1 2008
#"""
from ase import Atoms as ase_Atoms
from ase import Atom
from hotbit import Atoms
import numpy as nu
from numpy import pi
from math import sqrt, pi, atan, cos, sin
from box.mix import gcd
vec=nu.array
dot=nu.dot


    
def graphene(n1,n2,R,height=5.0):
    """
    Construct graphene lattice, multiply the primitive cell
    n1 x n2 times in corresponding directions.
        
         .-----.
        /     /
       /   X / a2  
      /     /
     X----->        
      a1            
    """
    a1=vec([R*nu.cos(pi/6)*2,0.,0.])
    a2=0.5*a1 + vec([0.,1.5*R,0.])
    #assert n2%2==0
        
    r=[]
    for i1 in range(n1):
        for i2 in range(n2):
            corner = i1*a1+i2*a2
            r.append(corner)
            r.append(corner+a1+vec([0.0,R,0.0]))
                
    cell=[[n1*a1[0], 0, 0],[n2*a2[0],n2*a2[1],0],[0,0,10]]                
    atoms=Atoms('C'*len(r),positions=r,cell=cell)
    atoms.center(vacuum=height/2,axis=2)
    atoms.set_pbc((True,True,False))
    return atoms

                         
def armchair_ribbon(n1,n2,R,pbc='z'):
    """
    Make ribbon out of graphene with armchair edges: n1 armchair periods in x-direction
    and n2 stacked armchairs.
    
    parameters:
    ===========
    n1:    length (units in periodic direction)
    n2:    width
    R:     bond distance
    pbc:   periodic direction ('z' or 'x')
    """
    a1=vec([3*R,0,0])
    a2=vec([0,2*R*nu.cos(pi/6),0])
    atoms=ase_Atoms()
    r=[]
    for i1 in range(n1):
        for i2 in range(n2):
            corner=i1*a1+i2*a2
            atoms += Atom('C',corner)
            atoms += Atom('C',corner+vec([R,0,0]))
            atoms += Atom('C',corner+vec([1.5*R,R*nu.cos(pi/6),0]))
            atoms += Atom('C',corner+vec([2.5*R,R*nu.cos(pi/6),0]))

    if pbc=='x':
        atoms.set_cell( [n1*3*R,n2*2*R*nu.cos(pi/6),1] )
        atoms.center(vacuum=5,axis=2)
        atoms.center(vacuum=5,axis=1)
        atoms.set_pbc((True,False,False))
        return atoms
    elif pbc=='z':
        atoms.translate( -atoms.get_center_of_mass() )
        atoms.rotate('z',nu.pi/2)
        atoms.rotate('x',nu.pi/2)
        atoms.center(vacuum=5,axis=1)
        atoms.translate( -atoms.get_center_of_mass() )
        zmin = atoms.get_positions()[:,2].min()
        atoms.translate( (0,0,-zmin) ) 
        atoms.set_cell( [n2*2*R*nu.cos(pi/6),1,n1*3*R] )
        atoms.set_pbc((False,False,True))
        return atoms
    else:
        raise NotImplementedError('pbc only along x or z')
            


    


def zigzag_ribbon(n1,n2,R,pbc='z'):
    """
    Make ribbon out of graphene with zigzag edges.
    """
    a1=vec([2*R*nu.cos(pi/6),0,0])
    a2=vec([0,3*R,0])
    r=[]
    for i1 in range(n1):
        for i2 in range(n2):
            corner=i1*a1+i2*a2
            r.append(corner+vec([R*nu.cos(pi/6),R/2,0]))
            r.append(corner+vec([0,R,0]))
            r.append(corner+vec([0,2*R,0]))
            r.append(corner+vec([R*nu.cos(pi/6),2.5*R,0]))
    cell=[n1*2*R*nu.cos(pi/6),n2*3*R,1]   
    elements=['C']*len(r)         
    atoms = ase_Atoms(elements,r,cell=cell)      
    if pbc=='x':
        atoms.set_cell( [n1*3*R,n2*2*R*nu.cos(pi/6),1] )
        atoms.center(vacuum=5,axis=2)
        atoms.center(vacuum=5,axis=1)
        atoms.set_pbc((True,False,False))
        return atoms
    elif pbc=='z':
        atoms.translate( -atoms.get_center_of_mass() )
        atoms.rotate('z',nu.pi/2)
        atoms.rotate('x',nu.pi/2)
        atoms.center(vacuum=5,axis=1)
        atoms.translate( -atoms.get_center_of_mass() )
        zmin = atoms.get_positions()[:,2].min()
        atoms.translate( (0,0,-zmin) ) 
        atoms.set_cell( [n2*2*R*nu.cos(pi/6),1,n1*3*R] )
        atoms.set_pbc((False,False,True))
        return atoms
    else:
        raise NotImplementedError('pbc only along x or z')
    


def nanotube_data(n,m,R=1.42):
    """
    Return a dictionary if miscellaneous nanotube data 
    """
    a=sqrt(3.0)*R
    data={}
    data['a'] = a
    data['diameter'] = a/pi*sqrt(n**2+m**2+n*m)
    data['radius'] = data['diameter']/2
    data['C'] = data['diameter']*pi
    d = gcd(n,m)
    if nu.mod(n-m,3*d)==0:
        dR = 3*d
    else:
        dR = d
    data['T'] = 3*R*sqrt(n**2+m**2+n*m)/dR
    data['metallic'] = (nu.mod(n-m,3)==0)
    data['semiconducting'] = not data['metallic']
    data['angle'] = atan( sqrt(3)*m/(m+2*n) )
    data['hexagons_per_cell'] = int(round(2*(m**2+n**2+m*n)/dR,0))
    if n==m:
        data['type'] = 'armchair'
    elif n==0 or m==0:
        data['type'] = 'zigzag'
    else:
        data['type'] = 'chiral'
    return data
    
    
def nanotube(n,m,R=1.42,chiral=True,element='C'):
    """
    Construct nanotube.
    
    @param n: chiral index n
    @param m: chiral index m
    @param R: nearest neighbor distance
    @param chiral: return hotbit.Atoms with container 'Chiral',
              using the minimum number of atoms in unit cell
    @param element: element symbol
        
    atoms.data will contain the dictionary for nanotube data.
    
    Reference: V. N. Popov, New J. Phys. 6, 17 (2004)
    """
    data = nanotube_data(n,m,R)
    T, radius, angle = data['T'], data['radius'], data['angle']
    if m>n:
        n, m = m, n
    L1, L2 = n, m
    dp = gcd(n,m)
    if nu.mod(L1-L2,3*dp)!=0:
        d = dp
    else:
        d = 3*dp
        
    N1 = (L1+2*L2)/d          # eq. (7)
    N2 = -(2*L1+L2)/d         # eq. (8)
    Nc = N1*L2-N2*L1          # eq. (13)
    phi1 = 2*pi*N2/Nc         # eq. (9)
    phi2 = -2*pi*N1/Nc        # eq. (10)
    t1 = T*L2/Nc              # eq. (11)
    t2 = -T*L1/Nc             # eq. (12)

    
    atoms = Atoms()
    lst = [(l,l) for l in range(m)] + [(l,m-1) for l in range(m,n)]
    for (l1,l2) in lst:
            phi = l1*phi1 + l2*phi2
            t   = l1*t1 + l2*t2
            # first atom is at origin, shifted by the symmetry operation in eq. (1)
            atoms += Atom('C',(radius*cos(phi),radius*sin(phi),t))
            # another basis atom is one-third on of operation S1*S2
            dphi = (phi1 + phi2)/3
            dt = (t1+t2)/3
            atoms += Atom('C',(radius*cos(phi+dphi),radius*sin(phi+dphi),t+dt))      
            
    atoms = Atoms(atoms,container='Chiral')
    angle, height = abs(phi2), abs(t2)
    atoms.set_container(angle=angle,height=height)
    if chiral:
        data['height']=height
        data['twist']=angle
    else:
        N = int( round(T/height) )
        cp = atoms.extended_copy((1,1,N))
        atoms = Atoms( cp,cell=(10,10,N*height),pbc=(False,False,True) )
        data['height']=N*height
        data['twist']=0.0
    atoms.data = data
    return atoms
    

def nanotube_old(el, a0, n, m, ncopy=1, verbose=False):
    '''
    Create a nanotube.
    
    @param el: element symbol
    @param a0: nearest neighbor distance
    @param n:  first index
    @param m:  second index
    @param ncopy: number of copies along z-axis
    @param verbose: print info related to carbon nanotubes 
    '''
    at = Atoms( pbc = ( False, False, True ) )

    sq3 = sqrt(3.0)

    gcn = gcd(n, m)
    
    a1 = nu.array( [ sq3/2,  0.5 ] ) * a0 * sq3
    a2 = nu.array( [ sq3/2, -0.5 ] ) * a0 * sq3

    h = float(float(n)-float(m))/float(3*gcn)

    if h-int(h) == 0.0:
        R = 3
    else:
        R = 1

    c = n*a1 + m*a2
    abs_c = sqrt(dot(c, c))

    a = ( -(2*m+n)*a1 + (2*n+m)*a2 )/(gcn*R)
    abs_a = sqrt(dot(a, a))

    eps = 0.01
    b = [ [ 1./3-eps, 1./3-eps ], [ 2./3-eps, 2./3-eps ] ]

    nxy = max(n, m)+100
    eps = 0.00001
    
    for x in xrange(-nxy, nxy):
        for y in xrange(-nxy, nxy):
            for b1, b2 in b:
                p = (x+b1)*a1 + (y+b2)*a2
                abs_p = sqrt(dot(p, p))

                sa = dot(p, a)/(abs_a**2)
                sc = dot(p, c)/(abs_c**2)

                if sa >= 0 and sa < 1-eps and sc >= 0 and sc < 1-eps:
                    r = ( cos(2*pi*sc)*abs_c/(2*pi), sin(2*pi*sc)*abs_c/(2*pi), sa*abs_a )
                    at += Atom( el, r ) 
    at.set_cell( ( 2*abs_c/(2*pi), 2*abs_c/(2*pi), ncopy*abs_a ) )
    b = at.copy()

    for i in range(ncopy-1):
        b.translate( ( 0.0, 0.0, abs_a ) )
        for j in b:
            at += j
    at.center(axis=2)
    at.set_pbc((False,False,True))
    if verbose:
        print 'D=',nu.sqrt(3)*a0/nu.pi * nu.sqrt( n**2+m**2+n*m )
    return at
            
##print help(make_zigzag_ribbon)
#try:
    #cmd=sys.argv[1]
#except:
    #funcs=dir()
    #for f in funcs:
        #if f.find('make')>=0: 
            #exec('print help(%s)' %f)
    #sys.exit()
    
#cmd='r,cell,elements='+cmd
#print cmd
#exec(cmd)

#xyz_output(elements,r,cell,file='atoms.xyz')

if __name__=='__main__':
    import ase
    c=armchair_ribbon(10,10,1.42)
    ase.view(c)

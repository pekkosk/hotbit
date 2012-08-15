#!/usr/bin/env python
"""
    Make various atomic structures into xyz-files.
    
    P. Koskinen 31.1 2008
"""
from ase import Atoms as ase_Atoms
from ase import Atom
import numpy as np
from math import sqrt, atan, cos, sin
from box.mix import gcd
vec=np.array
dot=np.dot
pi = np.pi


    
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
    from hotbit import Atoms
    
    if not isinstance(R,float): R=R[0]
    a1=vec([R*np.cos(pi/6)*2,0.,0.])
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

def ZGNR(n,units=1,pbc='z',R=1.42):
    """
    Make an n-zigzag graphene nanoribbon.
    
    parameters:
    ===========
    n:      ribbon width (atomic rows)
    length: ribbon length (unit cells)
    pbc:    periodic direction, 'x' or 'z'
    R:      bond length 
    """
    from hotbit import Atoms
    a0 = R*np.sqrt(3)
    atoms = ase_Atoms()
    for i in range(n):
        x = 1        
        if np.mod(i,2)==1:
            x=-1
        atoms += Atom('C',(0,3*R*i/2+3*R/4-x*R/4,0))
        atoms += Atom('C',(a0/2,3*R*i/2+3*R/4+x*R/4,0))
    
    a = atoms.copy()
    for i in range(units-1):
        b = a.copy()
        b.translate(((i+1)*a0,0,0))
        atoms += b
    if pbc=='x':
        atoms.set_pbc((True,False,False))
        atoms.set_cell((units*a0,1,1))
        atoms.center(vacuum=6,axis=1)
        atoms.center(vacuum=6,axis=2)
    elif pbc=='z':
        atoms.set_pbc((False,False,True))
        atoms.rotate('z',np.pi/2)
        atoms.rotate('x',np.pi/2)
        atoms.set_cell((n*3*R/2*1.2/2,1,units*a0))
        atoms.translate( -atoms.get_center_of_mass() )
        atoms.center(axis=2)
    return atoms
        
    

def AGNR(n,units=1,pbc='z',R=1.42):
    """
    Make an n-armchair graphene nanoribbon.
    
    parameters:
    ===========
    n:      ribbon width (atomic rows)
    length: ribbon length (unit cells)
    pbc:    periodic direction, 'x' or 'z'
    R:      bond length 
    """
    from hotbit import Atoms
    a = R*np.sqrt(3)
    atoms = ase_Atoms()
    for i in range(n):
        x0 = 0.0
        if np.mod(i,2)==1:
            x0 = 1.5*R
        atoms += Atom('C',(x0+0,i*a/2,0))
        atoms += Atom('C',(x0+R,i*a/2,0))
        
    a = atoms.copy()
    for i in range(units-1):
        b = a.copy()
        b.translate(((i+1)*3*R,0,0))
        atoms += b
    if pbc=='x':
        atoms.set_pbc((True,False,False))
        atoms.set_cell((units*3*R,1,1))
        atoms.center(vacuum=6,axis=1)
        atoms.center(vacuum=6,axis=2)
    elif pbc=='z':
        atoms.set_pbc((False,False,True))
        atoms.rotate('z',np.pi/2)
        atoms.rotate('x',np.pi/2)
        atoms.set_cell((1,1,units*3*R))
        atoms.translate( -atoms.get_center_of_mass() )
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
    from hotbit import Atoms
    a1=vec([3*R,0,0])
    a2=vec([0,2*R*np.cos(pi/6),0])
    atoms=ase_Atoms()
    r=[]
    for i1 in range(n1):
        for i2 in range(n2):
            corner=i1*a1+i2*a2
            atoms += Atom('C',corner)
            atoms += Atom('C',corner+vec([R,0,0]))
            atoms += Atom('C',corner+vec([1.5*R,R*np.cos(pi/6),0]))
            atoms += Atom('C',corner+vec([2.5*R,R*np.cos(pi/6),0]))

    if pbc=='x':
        atoms.set_cell( [n1*3*R,n2*2*R*np.cos(pi/6),1] )
        atoms.center(vacuum=5,axis=2)
        atoms.center(vacuum=5,axis=1)
        atoms.set_pbc((True,False,False))
        return atoms
    elif pbc=='z':
        atoms.translate( -atoms.get_center_of_mass() )
        atoms.rotate('z',np.pi/2)
        atoms.rotate('x',np.pi/2)
        atoms.center(vacuum=5,axis=1)
        atoms.translate( -atoms.get_center_of_mass() )
        zmin = atoms.get_positions()[:,2].min()
        atoms.translate( (0,0,-zmin) ) 
        atoms.set_cell( [n2*2*R*np.cos(pi/6),1,n1*3*R] )
        atoms.set_pbc((False,False,True))
        return atoms
    else:
        raise NotImplementedError('pbc only along x or z')
            


    


def zigzag_ribbon(n1,n2,R,pbc='z'):
    """
    Make ribbon out of graphene with zigzag edges.
    """
    from hotbit import Atoms
    a1=vec([2*R*np.cos(pi/6),0,0])
    a2=vec([0,3*R,0])
    r=[]
    for i1 in range(n1):
        for i2 in range(n2):
            corner=i1*a1+i2*a2
            r.append(corner+vec([R*np.cos(pi/6),R/2,0]))
            r.append(corner+vec([0,R,0]))
            r.append(corner+vec([0,2*R,0]))
            r.append(corner+vec([R*np.cos(pi/6),2.5*R,0]))
    cell=[n1*2*R*np.cos(pi/6),n2*3*R,1]   
    elements=['C']*len(r)         
    atoms = ase_Atoms(elements,r,cell=cell)      
    if pbc=='x':
        atoms.set_cell( [2*R*np.cos(pi/6),n1*3*R,1] )
        atoms.center(vacuum=5,axis=2)
        atoms.center(vacuum=5,axis=1)
        atoms.set_pbc((True,False,False))
        return atoms
    elif pbc=='z':
        atoms.translate( -atoms.get_center_of_mass() )
        atoms.rotate('z',np.pi/2)
        atoms.rotate('x',np.pi/2)
        atoms.center(vacuum=5,axis=1)
        atoms.translate( -atoms.get_center_of_mass() )
        zmin = atoms.get_positions()[:,2].min()
        atoms.translate( (0,0,-zmin) ) 
        atoms.set_cell( [n2*2*R*np.cos(pi/6),1,n1*3*R] )
        atoms.set_pbc((False,False,True))
        return atoms
    else:
        raise NotImplementedError('pbc only along x or z')
    


def nanotube_data(n,m,R=1.42):
    """
    Return a dictionary if miscellaneous nanotube data 
    """
    from hotbit import Atoms
    a=sqrt(3.0)*R
    data={}
    data['a'] = a
    data['diameter'] = a/pi*sqrt(n**2+m**2+n*m)
    data['radius'] = data['diameter']/2
    data['C'] = data['diameter']*pi
    if n==m:
        d = n
        dR = 3*n
    else:
        d = gcd(n,m)
        if np.mod(n-m,3*d)==0:
            dR = 3*d
        else:
            dR = d
        
    data['T'] = 3*R*sqrt(n**2+m**2+n*m)/dR
    data['metallic'] = (np.mod(n-m,3)==0)
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
    
    
    

def chiral_nanotube(n,m,R=1.42,element='C'):
    """
    Construct a nanotube with a chiral container.
    
    parameters:
    ===========
    n,m:     chiral indices
    R:       bond length
    element: element type 
    """
    from hotbit import Atoms
    a = np.sqrt(3)*R
    a1 = np.array([a,0,0])
    a2 = np.array([0.5*a,-1.5*R,0])
    
    rl = []
    shift = (a/2,-0.5*R,0)
    for i in range(n):
        origin = i*a1
        rl.append( origin )
        rl.append( origin+shift )
          
    for j in range(m):
        origin = (n-1)*a1 + (j+1)*a1
        rl.append( origin )
        rl.append( origin + shift )
        
    atoms = Atoms()
    C = n*a1 + m*a2
    Cn = np.linalg.norm(C)
    T = np.array([C[1],-C[0],0])
    t = T/np.linalg.norm(T)
    radius = Cn/(2*pi)
    
    atoms = Atoms()
    for r in rl:
        phi = np.dot(r,C)/Cn**2 * 2*pi
        atoms += Atom( element,(radius*np.cos(phi),radius*np.sin(phi),np.dot(r,t)) )
    
    atoms = Atoms(atoms,container='Chiral')
    height = np.abs( np.dot(a1-a2,t) )
    angle = -np.dot(a1-a2,C)/Cn**2 * 2*pi
    atoms.set_container(angle=angle,height=height) 
    
    data = nanotube_data(n,m,R)
    T, Ntot = data['T'],2*data['hexagons_per_cell']
    data['height']=height
    data['twist']=angle
    atoms.data = data
    return atoms
         


def nanotube(n,m,R=1.42,length=1,element='C'):
    '''
    Create a nanotube around z-axis.
    
    parameters:
    -----------
    n,m:    chiral indices
    R:      nearest neighbor distance
    length: number of unit cells
    element: element symbol
    '''
    from hotbit import Atoms
    at = Atoms( pbc = ( False, False, True ) )

    sq3 = sqrt(3.0)
    a0 = R
    gcn = gcd(n, m)
    
    a1 = np.array( [ sq3/2,  0.5 ] ) * a0 * sq3
    a2 = np.array( [ sq3/2, -0.5 ] ) * a0 * sq3

    h = float(float(n)-float(m))/float(3*gcn)

    if h-int(h) == 0.0:
        RR = 3
    else:
        RR = 1

    c = n*a1 + m*a2
    abs_c = sqrt(dot(c, c))

    a = ( -(2*m+n)*a1 + (2*n+m)*a2 )/(gcn*RR)
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
                    at += Atom( element, r ) 
    at.set_cell( ( 2*abs_c/(2*pi), 2*abs_c/(2*pi), length*abs_a ) )
    b = at.copy()

    for i in range(length-1):
        b.translate( ( 0.0, 0.0, abs_a ) )
        for j in b:
            at += j
    at.center(axis=2)
    rcm = at.get_center_of_mass()
    at.translate( (-rcm[0],-rcm[1],0) )
    at.set_pbc((False,False,True))
    at.data = nanotube_data(n,m)
    return at



def CNT_list(i=None):
    """
    Select a carbon nanotube from a list.
    
    Order all CNTs according to the number of atoms in the
    unit cell, with diameters larger than 2.0 Angstrom, 
    assuming 1.42 Angstrom bond length.
    
    
    parameters:
    -----------
    i:     return the (n,m) tuple of i'th CNT in the list.
           If i==None, return the whole list of 
           (N,d,L,n,m) -tuples (N=atom count, d=diameter,
           L=unit cell length, n,m=chiral indices)
    """
    R = 1.42
    dmin = 2.0
    
    data  = []
    for n in range(1,20):
        for m in range(n+1):
            cnt = nanotube_data(n,m,R=R)
            N = cnt['hexagons_per_cell']*2
            d = cnt['diameter']
            L = cnt['T']
            if d<dmin: break
            data.append( (N,d,L,n,m) )    
            
        #view( nanotube(n,m) )
        #raise SystemExit

    def comp(x,y):
        if x[0]==y[0]:
            if x[1]==y[1]:
                return cmp(x[2],y[2])
            else:
                return cmp(x[1],y[1])
        else:
            return cmp(x[0],y[0])

    data.sort(comp)
    if i==None:
        return data
    else:
        return (data[i][3],data[i][4])


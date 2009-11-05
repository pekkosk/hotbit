#"""
    #Make various atomic structures into xyz-files.
    
    #P. Koskinen 31.1 2008
#"""
from ase import Atoms,Atom
import numpy as nu
from numpy import pi
vec=nu.array
sqrt=nu.sqrt
dot=nu.dot
cos=nu.cos
sin=nu.sin

#def xyz_output(elements,coords,cell,file='out.xyz'):
    #f=open(file,'w')
    #dr=ranges(coords)
    #print>>f, len(elements)
    ##print>>f, 'UnitCell: ',mix.a2s(cell[0]),mix.a2s(cell[1]),mix.a2s(cell[2]),'dr:',dr[0],dr[1],dr[2]
    ##print>>f, 'UnitCell: ',cell[0][0],cell[1][1],cell[2][2],'90 90 90'
    #print>>f, 'UnitCell: ',cell[0],cell[1],cell[2],'90 90 90'
    #for e,r in zip(elements,coords):
        #print>>f, e,r[0],r[1],r[2]
    #f.close()
    
#def ranges(coords):
    #""" Return the x,y and z ranges of the coordinates. dx,dy,dz."""
    #mn=coords.min(axis=0)
    #mx=coords.max(axis=0)
    #return mx-mn
    
#def not_implemented():
    #raise NotImplementedError('todo...')
    
#def make_fcc(nx,ny,nz,element,a=None,R=None):
    #"""
    #"""
    #if a==R==None or (a!=None and R!=None):
        #raise ValueError('either a or R has to be specified; not both')
    #if a==None:
        #a=nu.sqrt(2.0)*R
    #else:
        #R=a/nu.sqrt(2.0)
    #ax=vec([a,0,0])
    #ay=vec([0,a,0])
    #az=vec([0,0,a])
    #r=[]
    #for ix in range(nx):
        #for iy in range(ny):
            #for iz in range(nz):       
                #corner=ix*ax+iy*ay+iz*az
                #r.append(corner)
                #r.append(corner+0.5*(ax+ay))
                #r.append(corner+0.5*(ax+az))
                #r.append(corner+0.5*(ay+az))
    #cell=[nx*ax,ny*ay,nz*az]
    #elements=[element]*len(r)
    #return vec(r),cell,elements
                
    
#def make_fcc_tube():
    #not_implemented()

#def make_bcc():
    #not_implemented()
    
#def make_sc():
    #not_implemented()
    
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
                         
def armchair_ribbon(n1,n2,R):
    """
    Make ribbon out of graphene with armchair edges: n1 armchair periods in x-direction
    and n2 stacked armchairs.
    """
    a1=vec([3*R,0,0])
    a2=vec([0,2*R*nu.cos(pi/6),0])
    atoms=Atoms()
    r=[]
    for i1 in range(n1):
        for i2 in range(n2):
            corner=i1*a1+i2*a2
            atoms += Atom('C',corner)
            atoms += Atom('C',corner+vec([R,0,0]))
            atoms += Atom('C',corner+vec([1.5*R,R*nu.cos(pi/6),0]))
            atoms += Atom('C',corner+vec([2.5*R,R*nu.cos(pi/6),0]))
            #r.append(corner)
            #r.append(corner+vec([R,0,0]))
            #r.append(corner+vec([1.5*R,R*nu.cos(pi/6),0]))
            #r.append(corner+vec([2.5*R,R*nu.cos(pi/6),0]))
    atoms.set_cell( (n1*3*R,n2*2*R*nu.cos(pi/6),0) )
    atoms.center(vacuum=2,axis=2)
    return atoms


def zigzag_ribbon(n1,n2,R):
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
    cell=[n1*2*R*nu.cos(pi/6),n2*3*R,0]   
    elements=['C']*len(r)         
    return Atoms(elements,r,cell=cell)



def gcd(a, b):
    if b == 0:
        return a
    if a == 0:
        return b
    
    h = min(a, b)
    while float(a)/float(h)-int(a/h) != 0 or float(b)/float(h)-int(b/h) != 0:
        h -= 1

    return h


def nanotube(el, a0, n, m, ncopy=1, verbose=False):
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

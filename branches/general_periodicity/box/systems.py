#"""
    #Make various atomic structures into xyz-files.
    
    #P. Koskinen 31.1 2008
#"""
from ase import Atoms
import numpy as nu
from numpy import pi
vec=nu.array


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
    
def graphene(n1,n2,R):
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
    assert n2%2==0
        
    r=[]
    for i1 in range(n1):
        for i2 in range(n2):
            corner = i1*a1+i2*a2
            r.append(corner)
            r.append(corner+a1+vec([0.0,R,0.0]))
                
    cell=( n1*a1[0], n2*a2[1], 0 )                
    atoms=Atoms('C'*len(r),positions=r,cell=cell)
    atoms.center(vacuum=5,axis=2)
    atoms.set_pbc((True,True,False))
    return atoms
                         
def armchair_ribbon(n1,n2,R):
    """
    Make ribbon out of graphene with armchair edges: n1 armchair periods in x-direction
    and n2 stacked armchairs.
    """
    a1=vec([3*R,0,0])
    a2=vec([0,2*R*nu.cos(pi/6),0])
    r=[]
    for i1 in range(n1):
        for i2 in range(n2):
            corner=i1*a1+i2*a2
            r.append(corner)
            r.append(corner+vec([R,0,0]))
            r.append(corner+vec([1.5*R,R*nu.cos(pi/6),0]))
            r.append(corner+vec([2.5*R,R*nu.cos(pi/6),0]))
    cell=[n1*3*R,n2*2*R*nu.cos(pi/6),0]            
    elements=['C']*len(r)
    return vec(r),cell,elements


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

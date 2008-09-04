from hotbit.parametrization import KSAllElectron
import pylab as pl
import numpy as nu
from box import mix
from box.interpolation import Function

integrals =['dds','ddp','ddd','pds','pdp','pps','ppp','sds','sps','sss']

def plot_table(parfile,screen=False,s1=None,s2=None,der=0):
    """ Plot table. """
    if s1==None or s2==None:
        s1,s2=parfile.split('.')[0].split('_')        
    if s1==s2:
        nel=1
        pairs=[(s1,s2)]
    else:
        nel=2
        pairs=[(s1,s2),(s2,s1)]        
    
    pl.rc('figure.subplot',wspace=0.0001)
    pl.rc('figure.subplot',hspace=0.0001)
    rgrid, tables=read_table(parfile,s1,s2)
    mx=max( tables[0].flatten() )
    
    for i in range(10):
        name=integrals[i]
        ax=pl.subplot(5,2,i+1)
        for p,(e1,e2) in enumerate(pairs):
            if p==1: s='--'
            else: s='-'
            if der==0:
                grid=rgrid
                H=tables[p][:,i]
                S=tables[p][:,i+10]
            elif der==1:
                grid=nu.linspace(rgrid[0],rgrid[-1],3*len(rgrid))
                hf=Function('spline',rgrid,tables[p][:,i])
                sf=Function('spline',rgrid,tables[p][:,i+10])
                H=nu.array([hf(r,der=1) for r in grid])
                S=nu.array([sf(r,der=1) for r in grid])
            pl.plot(grid,H,c='r',ls=s,label='%s%s: H' %(e1,e2))
            pl.plot(grid,S,c='b',ls=s,label='%s%s: S' %(e1,e2))                
            pl.axhline(c='k',ls='--')
            pl.title(name,position=(0.9,0.8)) 
            if ax.is_last_row():
                pl.xlabel('r (Bohr)')                                        
            else:
                pl.xticks([],[])
            if not ax.is_first_col():                   
                pl.yticks([],[])
            pl.ylim(-mx,mx)
            pl.xlim(0)
            pl.legend(loc='upper left')
    if screen:                
        pl.show()
    else:
        name='%s_%s_par.png' %(s1,s2)
        if der==1:
            name='der_'+name
        pl.savefig(name)   
        
        
def compare_tables(parfile1,parfile2,s1=None,s2=None,screen=False):
    """ Plot table. """
    if s1==None or s2==None:
        s1,s2=parfile1.split('.')[0].split('_')        
    if s1==s2:
        nel=1
        pairs=[(s1,s2)]
    else:
        nel=2
        pairs=[(s1,s2),(s2,s1)]        
    
    pl.rc('figure.subplot',wspace=0.0001)
    pl.rc('figure.subplot',hspace=0.0001)
    rgrid1, tables1=read_table(parfile1,s1,s2)
    rgrid2, tables2=read_table(parfile2,s1,s2)
    mx=max( tables1[0].flatten() )*0.5
    
    for i in range(10):
        name=integrals[i]
        ax=pl.subplot(5,2,i+1)
        for p,(e1,e2) in enumerate(pairs):
            if p==1: s='--'
            else: s='-'
            # first table
            pl.plot(rgrid1,tables1[p][:,i],lw=5,c='r',alpha=0.3,ls=s,label='%s%s: H' %(e1,e2))
            pl.plot(rgrid1,tables1[p][:,i+10],lw=5,alpha=0.3,c='b',ls=s,label='%s%s: S' %(e1,e2))
            
            # second table
            pl.plot(rgrid2,tables2[p][:,i],lw=2,c='r',ls=s,label='%s%s: H' %(e1,e2))
            pl.plot(rgrid2,tables2[p][:,i+10],lw=2,c='b',ls=s,label='%s%s: S' %(e1,e2))
            pl.axhline(c='k',ls='--')
            pl.title(name,position=(0.9,0.8)) 
            if ax.is_last_row():
                pl.xlabel('r (Bohr)')                                        
            else:
                pl.xticks([],[])
            if not ax.is_first_col():                   
                pl.yticks([],[])
            pl.ylim(-mx,mx)
            pl.xlim(0)
            #pl.legend(loc='upper left')
    if screen:                
        pl.show()
    else:
        pl.savefig('%s_%s_comparison.png' %(s1,s2))           
        pl.close()
        
        
def read_table(parfile,s1,s2):
    """ Read parameter table from file parfile for elements with symbols s1 and s2. 
    
    return list of tables [s1_s2_table,s2_s1_table] (or only other if s1==s2)
    """
    f=open(parfile)
    nel=[1,2][s1==s2]
    tab=mix.find_value(parfile,'%s_%s_table' %(s1,s2),fmt='matrix')
    rgrid=tab[:,0]
    table=[tab[:,1:]]
    if s1!=s2:
        tab=mix.find_value(parfile,'%s_%s_table' %(s2,s1),fmt='matrix')
        table.append(tab[:,1:])
    f.close()                
    return rgrid, table
    
    
def tail_smoothening(x,y):
    """ For given grid-function y(x), make smooth tail.
    
    Aim is to get (e.g. for Slater-Koster tables and repulsions) smoothly
    behaving energies and forces near cutoff region.
    
    Make is such that y and y' go smoothly exactly to zero at last point.
    Method: take largest neighboring points y_k and y_(k+1) (k<N-3) such
    that line through them passes zero below x_(N-1). Then fit
    third-order polynomial through points y_k, y_k+1 and y_N-1.
    
    Return:
    smoothed y-function on same grid.
    """
    if all(abs(y)<1E-10):
        return y
    N=len(x)
    xmax=x[-1]
    for i in range(N-3,1,-1):
        x0i=x[i]-y[i]/( (y[i+1]-y[i])/(x[i+1]-x[i]) )
        if x0i<xmax:
            k=i
            break
    if k<N/4:
        for i in range(N):
            print x[i],y[i]
        raise RuntimeError('Problem with tail smoothening: requires too large tail.')                
    if k==N-3:
        y[-1]=0.0
        return y
    else:        
        # g(x)=c2*(xmax-x)**m + c3*(xmax-x)**(m+1) goes through (xk,yk),(xk+1,yk+1) and (xmax,0)
        # Try different m if g(x) should change sign (this we do not want)
        sgn=nu.sign(y[k])
        for m in range(2,10):
            a1, a2=(xmax-x[k])**m, (xmax-x[k])**(m+1)
            b1, b2=(xmax-x[k+1])**m, (xmax-x[k+1])**(m+1)
            c3=(y[k]-a1*y[k+1]/b1)/(a2-a1*b2/b1)
            c2=(y[k]-a2*c3)/a1
            for i in range(k+2,N):
                y[i]=c2*(xmax-x[i])**2 + c3*(xmax-x[i])**3
            y[-1]=0.0 #once more excplicitly            
            if all(y[k:]*sgn>=0): 
                break
            if m==9:
                raise RuntimeError('Problems with function smoothening; need for new algorithm?')
    return y            
            
def IP_EA(symb,remove_orb,add_orb,remove,add):
    """ Return ionization potential and electron affinity for given atom. 
    
    parameters: 
    -----------
    symb: element symbol
    remove_orb: orbital from where to remove atoms (e.g. '2p')
    add_orb: orbital from where to add atoms (e.g. '2p')
    remove: how many electrons to remove
    add: how many electrons to add
         (remove and add can be different from 1.0 if DFT should not
         be stable to e.g. adding one full electron)
         
    Fit second order curve for 3 points and return IP and EA for full
    electron adding and removal.             
    """
    from box.data import atom_occupations
    # neutral atom
    atom=KSAllElectron(symb)
    atom.run()
    e0=atom.get_energy()
    # remove electrons -> positive ion
    occu=atom_occupations[symb].copy()
    occu[remove_orb]-=remove
    ep=KSAllElectron(symb,occu=occu)
    ep.run()
    ep=ep.get_energy()-e0
    # add electrons -> negative ion
    occu=atom_occupations[symb].copy()
    occu[add_orb]+=add
    en=KSAllElectron(symb,occu=occu)
    en.run()
    en=en.get_energy()-e0
    # e(x)=e0+c1*x+c2*x**2 =energy as a function of additional electrons
    c2=(en+ep*add/remove)/(add*(remove+add))
    c1=(c2*remove**2-ep)/remove
    IP=-c1+c2
    EA=-(c1+c2)
    return IP, EA
    
    
    
    
    return (e1-e0)/electrons        
        
def ionization_potential(symb,remove,electrons=1.0):
    """ Return ionization potential of given atom. 
    
    parameters: 
    -----------
    symb: element symbol
    remove: orbital from where electron is removed (e.g. '2p')
    electrons: how many electrons to remove. Can be fractional number
               if DFT should not be stable. IP is scaled by
               electrons^-1 in the end to extrapolate to electrons=1.
    
    """
    from box.data import atom_occupations
    occu=atom_occupations[symb].copy()
    # neutral atom
    atom=KSAllElectron(symb)
    atom.run()
    e0=atom.get_energy()
    # negative ion
    occu[remove]-=electrons
    ion=KSAllElectron(symb,occu=occu)
    ion.run()
    e1=ion.get_energy()
    return (e1-e0)/electrons

def electron_affinity(symb,add,electrons=1.0):
    """ Return electron affinity of given atom. 
    
    parameters: 
    -----------
    symb: element symbol
    add: orbital where electron is added (e.g. '2p')
    electrons: how many electrons are added. Can be fractional number
               if DFT should not be stable. EA is scaled by
               electrons^-1 in the end to extrapolate to electrons=1.
    """
    from box.data import atom_occupations
    occu=atom_occupations[symb].copy()
    # neutral atom
    atom=KSAllElectron(symb)
    atom.run()
    e0=atom.get_energy()
    # positive ion
    occu[add]+=electrons
    ion=KSAllElectron(symb,occu=occu)
    ion.run()
    e1=ion.get_energy()
    return (e0-e1)/electrons
    
    
                
        
if __name__=='__main__':
    #plot_table('Au_Au.par',s1='Au',s2='Au',screen=True,der=0)
    #compare_tables('Au_Au.par','Au_Au_NR.par',s1='Au',s2='Au',screen=False)                
    
    x=nu.array([1,2,3,4,5,6,7,8,9,10])
    y=nu.array([100,70,30,10,5,2,0.5,0.1,0.05,0.0001])
    pl.plot(x,y)
    y=tail_smoothening(x,y)
    pl.plot(x,y)
    pl.show()
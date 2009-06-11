#!/usr/bin/env python

# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.
'''
    A module containing miscellanous utility functions.

    Author: P. Koskinen 11.9 2006-


'''
import numpy as nu
import time

def kronecker(i,j):
    if i==j:
        return 1
    else:
        return 0
    
def matrix_pprint(m,fmt='%.4f'):
    """ pretty-pring matrix. """
    assert len(m.shape)==2
    for i in range(len(m[:,0])):
        print a2s(m[i,:],fmt=fmt)

class Timer:
    def __init__(self):
        self.t0=time.time()
        self.t=time.time()
        
    def __call__(self):
        t=time.time()
        dt=t-self.t
        self.t=t
        print 'dt',dt     
        
    def __del__(self):
        t=time.time()
        print 'total time',t-self.t0        
    

def print_timing(func):
    def wrapper(*arg):
        t1 = time.time()
        res = func(*arg)
        t2 = time.time()
        print '%s took %0.3f ms' % (func.func_name, (t2-t1)*1000.0)
        return res
    return wrapper

def random_direction(vector=False):
    """ Return random direction from 4*pi (phi,theta) or unit vector. """
    phi=nu.random.random()*2*nu.pi
    theta=nu.arccos(1-2*nu.random.random())
    if not vector:
        return phi,theta
    else:
        return nu.array([nu.sin(theta)*nu.cos(phi),nu.sin(theta)*nu.sin(phi),nu.cos(theta)])
    
def spherical_to_cartesian(theta,phi,r=1.0):
    """ Transform spherical coordinates to cartesian ones. """
    return nu.array([r*nu.sin(theta)*nu.cos(phi),r*nu.sin(theta)*nu.sin(phi),r*nu.cos(theta)])
    
def cosd(x): 
    return cos(x*2*pi/360.0)

def sind(x): 
    return sin(x*2*pi/360.0)

def base_and_extension(filename):
    """ Separate the basename and extension from given file name. """
    i=filename.rfind('.')
    return filename[:i],filename[i+1:]

def sec_to_time(secs):
    """return days,hours,minutes and seconds in list [dd,hh,mm,ss]"""
    d = int(secs // 86400)
    s = secs % 86400
    h = int(s // 3600)
    s = secs % 3600
    m = int(s // 60)
    s = secs % 60
    s = int(s+.5)
    return [d,h,m,s]

def clean_indentation(st):
    """ 
    Return multiple-line string into lowest indentation equal to zero. 
    Remove initial and trailing empty line.
    """
    lines=st.splitlines()
    mn=100
    for line in lines:
        ns=0
        if len(line)<=1:
            continue
        else:
            for c in line:
                if c==" ": ns+=1
                else: break
        mn=min(ns,mn)
    new=''
    for line in lines:
        new+=line[mn:]+'\n' 
    first=False
    for c in new:
        if c=='\n' and not first:
            new=new[1:]
        else:
            first=True    
    return new[:-1]

def type_check(x,t):
    from numpy import array
    tx=type(x)
    types={'float':type(1.0),'str':type(''),'list':type([]),'vector':type(array([])),'int':type(1)}
    if types[t]!=tx:
        error_exit('variable %s is not of type %s' %(str(x),str(tx)) )
    
def todo(s):
    """ Give todo-message. """
    mix.error_exit('The function %s is not implemented yet.')
    
def max_norm(vecs):
    """ Return the maximum norm of given list of vectors. """
    mx=0.0
    assert type(vecs)==type([])
    for vec in vecs:
        mx=max( mx,norm(vec) )
    return mx

def select_vector_mode(vec,mode):
    """ Change given vector into given mode (default or global)
    If mode is the same as vector, do nothing.
    """
    from numpy import array as vector
    from numpy import mod
    listvec=( type(vec)==type([]) )
    if listvec and mode=='default':
        return vec
    elif not listvec and mode=='global':
        return vec
    elif listvec and mode=='global':
        x=list(vec[0])
        for vec in vec[1:]:
            x=x+list(vec)
        return vector(x)        
    elif not listvec and mode=='default':
        if mod(len(vec),3)!=0: error_exit('Vector dimension not divisible by 3!')
        N=len(vec)/3
        x=[]
        for i in range(N):
            x.append( vector(vec[i*3:i*3+3]) )
        return x    
    

def eof(file,tol=3):
    """ Test enf of file within given byte tolerance. """
    from os import lseek
    now=file.tell()
    end=lseek(file.fileno(),0,2)
    file.seek(now)
    if end-now<=tol:
        return True
    else:
        return False  
    
def a2s(a,fmt='%.14g',end=''):
    """ Transforms vector array into a string with given format."""
    st=''
    for i in range( len(a) ):
        st=st+str(fmt%a[i]+'\t'  )
    return st+end
    
def execute(cmd,echo=True):
    from os import system
    from sys import exit
    if echo: 
        s=cmd
        if len(s)>80:
            s=s.split()
            s1=' '.join(s[:-1])
            print 'Executing:',s1,'...'
            print '       ...',s[-1],'...'
        else: print 'Executing:',s,'...'
    system(cmd)
    #system(d    
    #f=popen(cmd)
    #fail=f.readline()
    #if fail.strip()=='' or fail.strip()=='0':
        #pass
    #else:
        #if echo: print 'Executing "',cmd,'" failed. Exit now.'
        #exit(1)

def error_exit(s):
    """ Simpy exit routine. """
    from sys import exit
    print '\n \n Error:'
    #s=s.splitlines()
    #for l in s: 
        #print '    ',l.strip()
    #print '\n\n'
    raise RuntimeError(s)
    #exit(1)

def norm(a):
    """
    Return the norm of the vector a.
    """
    from numpy import sqrt,dot
    return sqrt( dot(a,a) )

def gauss_fct(x,mean=0.,sigma=1.):
    """
    Returns 1/sqrt(2*pi*sigma**2)*exp(-(x-mean)**2/(2*sigma**2)
    mean=0 and sigma=1 on default
    """
    from numpy import sqrt,pi,exp
    return 1./sqrt(2*pi*sigma**2)*exp(-(x-mean)**2/(2*sigma**2) )

def lorenzian(x,mean,width):
    """ Return normalized Lorenzian with given mean and broadening. """
    return (width/nu.pi)/((x-mean)**2+width**2)

def broaden(x,y=None,width=0.05,function='gaussian',extend=False,N=200,a=None,b=None):
    """ 
    Broaden a peaked distribution (DOS,optical spectra,...).
    
    parameters:
    -----------
    x: data points (~energy axis)
    y: heights of the peaks given with x. Default is one for all.
    width: width parameter specific for given broadening function.
    function: 'gaussian' or 'lorenzian'
    extend: if True, extend xrange bit beyond min(x) and max(x) (unless [a,b] given)
    N: number of points in output
    a: if defined, is used as the lower limit for output
    b: if defined, is used as the upper limit for output
    
    return: xgrid, broadened distribution
    """
    if function=='lorenzian':
        raise NotImplementedError('Only Gaussian works now for efficiency reasons')
    
    if y==None:
        y=nu.ones_like(x)
    dx=[0.0,4*width][extend]
    if a==None:
        mn=min(x)-dx
    else:
        mn=a
    if b==None:
        mx=max(x)+dx
    else:
        mx=b
        
    xgrid = nu.linspace(mn,mx,N)
    ybroad= nu.zeros_like(xgrid)
    for xi,yi in zip(x,y):
        h = yi*nu.exp( -(xgrid-xi)**2/(2*width**2) ) / (nu.sqrt(2*nu.pi)*width)
        ybroad = ybroad + h
    return xgrid, ybroad 
    

def grid(min,max,N):
    """ 
    Returns a grid with min and max as end-points and (N-1) divisions.
    """
    from numpy import arange
    return min+arange(N)/(N-1.0)*(max-min)

def true_false(s):
    """
    Interpret string s as boolean ans return True or False.
    """
    s2=s.lower().strip()
    if( s2=='true' or s2=='t' or s2=='yes' or s2=='y' ):
        return True
    elif( s2=='false' or s2=='f' or s2=='no' or s2=='n' ):
        return False
    else:
        error_exit('String "'+s2+'" could not be interpreted as boolean')
        
    
def find_value(inp,key,fmt='default',default=None,position='start'):
    '''
    Find value for key from file.
    
    Usage: value_string=find_key(f,key)
           value_string_list=find_key(f,key,fmt='default')
    
    Parameters:
    -----------
      inp - file object or file name
      key  - key string (corresponding to 'key=value')
      fmt  - the format of the value:
             'default' = return first value as string from the same line
             'all'     = return the whole line from the same line
             'onestr'  = return the whole line from the same line as one string
             'matrix'  = return a matrix coming after key as float
             'strings' = return the non-empty lines coming after the key as str
             'bool'    = interpret the line as boolean
             'test'    = return True if key found, otherwise False
       default - returned value if key not found
       position - position to start searching
             'start'   = beginning of the file object
             'current' = current position
    Example: f=open('data.txt')
             mass=float( find_value(f,'mass') )
             mat = find_value(f,'m','matrix')
    '''
    f,opened=file_safeopen(inp)
    if position=='start': 
        f.seek(0)
    else:
        pass       
    empty=0
    ret=None
    while 1:
        line=f.readline()
        if(len(line)==0): empty+=1
        if(empty>1000): break
        i=line.find('=')
        hlp=line.split('=')
        if hlp[0].strip()==key.strip() and i>=0: 
            if fmt=='test':    ret=True 
            if fmt=='bool':    ret=true_false(hlp[1])
            if fmt=='matrix':  ret=read(f) 
            if fmt=='strings': ret=read(f,fmt='string')
            if len(hlp)>1: value=hlp[1].split() 
            if fmt=='default': ret=value[0]
            if fmt=='onestr':  ret=hlp[1][:-1]
            elif fmt=='all':   ret=value[0:]
            break
    if opened: f.close()
    if fmt=='test' and ret==None:
        return False
    elif ret!=None:
        if fmt=='strings' and type(ret)!=type([]):
            ret=[ret]
        return ret
    elif default!=None:
        return default
    else:
        raise RuntimeError('key '+key+' not found from file '+f.name)


def read_file(file,mode='all',commentchar='#'):
    """
    Read the lines in file into list of strings.
    mode:   all         - read the whole file
            comments    - read only comments
            nocomments  - everything but comments
    """
    from os.path import isfile

    opened=False
    if isinstance(file,str):
        if isfile(file):
            f=open(file,'r')
            opened=True
        else: error_exit('\n File '+file+' does not exist!')
    else:
        f=file
        f.seek(0)
    lines=f.readlines()
    if mode=='all':
        return lines
    else:
        ret=[]
        for line in lines:
            line=line.lstrip()
            if len(line)==0: continue
            if line[0]==commentchar and mode=='comments':   ret.append(line)
            if line[0]!=commentchar and mode=='nocomments': ret.append(line)
        return ret
    if opened: f.close()   
     
def file_safeopen(file,action='r'):
    """
    Return file if file is a file-object, otherwise return the opened file 
    object and True/False if the file was opened or not. Action 'r','w','rw'.
    """
    from os.path import isfile
    opened=False
    if isinstance(file,str):
        if isfile(file):
            f=open(file,action)
            opened=True
        else: 
            error_exit('\n File '+file+' does not exist!')
    else:
        f=file
    return (f,opened)
    
def identify_column(col,file,commentchar='#'):
    """
    Return the column index (starting from 0), identified by col.
    The first row identifies the columns, and starts with commentchar.
    Example: i=mix.identify_column('energy','stuff.out') -> i=2 with
    stuff.out begining e.g. # time temperature energy distance ...
    """
    from string import maketrans
    f,opened = file_safeopen(file)
    f.seek(0)
    columns=f.readline()
    i = columns.find(col)
    if i<0:
        error_exit('[identify_column] column '+col+' not found in '+file)
    columns=columns.translate( maketrans('#.,;:|','      ') ).split()
    if opened: f.close()

    return columns.index(col)


def read_column(col,file,commentchar='#'):
    """
    Load a column with name col into one-dimensional NumPy array.
    The first row identifies the columns, and starts with commentchar.
    Example: energies=mix.read_column('energy','stuff.out')
    (stuff.out begins e.g. # time temperature energy distance ...)
    """
    from numpy import array
    import sys
    f,opened = file_safeopen(file)
    dat = read(f)
    if opened: f.close()
    return dat[:,identify_column(col,file)]
    
    
def read(file,commentchar='#',fmt='array'):
    """ 
    Read a data set from given file.
    
    * file -- file name of file object (for object the "next" set is read)
    * commentchar -- character for comments, which are not read
    * fmt -- format for output:
        
        - 'array' -- return two-dimensional numpy array
        - 'string' -- return the set as a list of strings
        - 'lines' -- return following non-empty lines as one string
    """
    from os.path import isfile
    if fmt=='array': from numpy import array
    f,opened = file_safeopen(file)
    r=[]
    col=-1
    start=True
    while 1:        
        line = f.readline()
        if not line: break              #end of file
        line = line.lstrip().strip()
        if line.isspace() or len(line)==0:
            if start==True: continue    # initial blank line
            elif start==False: break    # next blank line stops the data set
        if line[0]=="#": continue       # ignore comment line
        cl   = len(line.split())
        #if col!=-1 and cl!=col: 
            #break # different amount of columns->stop
        if fmt=='array':
            r.append([float(x) for x in line.split()])
        elif fmt=='string':
            if line[-1]=='\n': line=line[:-1]
            r.append(line)       
        else: error_exit('[mix.read] Invalid format.')
        col = cl
        start=False
    if opened: f.close()
    if r==[]: return None
    if fmt=='array': 
        return array(r)
    elif fmt=='lines':
        return '\n'.join(r)
    else: 
        return r
    

def write(file,a,fmt='%g'):
    """
    Write a two-dimensional NumPy array a in a tabular form.
    file can be file name or file object.
    fmt is the format of output (Default: fmt='%g')
    Example:
      write('data.out',array,'%14.2f')
      write(f,r)
    
    """
    try:
        f=open(file,'w')
    except:
        f=file
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            f.write(fmt %a[i,j]+'\t')
        f.write('\n')
    f.flush()
               
def abs_sum(a):
    """ Return the sum of all elements (absolute values) in array. """
    return abs(a).sum()
      
      
def gplot(string,file=None,delete=False,view=False):
    """Make a gnuplot file and call gnuplot to make an image from string.
       If delete=True, remove the intermediate .gp -files
    """
    from os import system
    from os import environ
    
    if file==None:
        for line in string.splitlines():
            if line.find('set output')>=0: 
                psfile=line.split("'")[1]
        file=psfile.split('.')[0]+'.gp'
    f=open(file,'w')
    for line in string.splitlines():
        f.write(line.lstrip()+'\n')
    f.close()
    execute('gnuplot %s' %file, echo=False)
    if delete: system('rm %s' %file)
    psviewer=environ.get('PSVIEWER')
    if psviewer=='': psviewer='gs'
    i = file.rfind('.')
    if i<0:
        psfile=file+'.ps'
    else:
        psfile=file[:i]+'.ps'
    if view: execute('%s %s&' %(psviewer,psfile) )
    
    
def add_to_file(s,file):
    """Add string to a file."""
    f=open(file,'r')
    lines=f.readlines()
    f.close()
    f=open(file,'w')
    lines.append(s+'\n')
    f.writelines(lines)
    f.close()


def add_file_appendix(file,app):
    """
    Adds appendix to file name, preserving the file type. 
    Examples: 
        'result.out' + '1' --> 'result_1.out'
        'result' '1' --> 'result_1'
        'result.out.xyz' '1' --> 'result.out_1.xyz'        
    """
    i = file.rfind('.')
    if i<0:
        return file+'_'+app
    else:
        return file[:i]+'_'+app+file[i:]
 #-- SplineFunction(x,y) (y 1-or multidim) f(x,der=0)


def parse_name_for_atoms(atoms):
    """
    Returns a name for the atoms based on the symbols.
    """
    symbols = atoms.get_chemical_symbols()
    dict = {}
    for symbol in symbols:
        n = symbols.count(symbol)
        if n not in dict:
            dict.update({symbol:n})
    name = ''
    for symbol in dict:
        if dict[symbol] < 2:
            name += symbol
        else:
            name += symbol+str(dict[symbol])
    return name


class AnalyticFunction:
    """ Class for defining analytic functions with strings."""
    
    def __init__(self,expr,variable='x',**kwargs):
        """ e.g. f=AnalyticFunction('a*x**2+b*x',a=1,b=2). free=independent variable """
        self.f=expr
        self.args=kwargs
        self.args['variable']=variable
        
    def __call__(self,variable,**kwargs):
        """ variable is the independent variable """
        self.args.update(kwargs)
        self.args.update({self.args['variable']:variable})
        f=eval(self.f,globals(),self.args)
        return f
   
    def set(self,**kwagrs):
        self.args.update(kwargs)
        
    def info(self):
        print 'function:',self.f
        print 'current arguments:',self.args
        
    def plot(self,a=None,b=None,N=200,out='screen'):
        import pylab as pl
        if a==b==None:
            a,b=(-10,10)
        
        x=nu.linspace(a,b,N)
        f=[self(tx) for tx in x]
        pl.plot(x,f)
        if out=='screen':
            pl.show()
        else:
            pl.savefig(out)
                    
    
class IFunction:
    """Simple class for interpolating functions on a grid."""    
    def __init__(self,x,y,s=0):
        """ x is the grid and y the corresponding values.
            s is the "smoothness", default s=0 (f(x_i)=y_i exactly)
        """        
        from scipy.interpolate import splrep 
        from scipy.interpolate import splev
        from scipy.interpolate import splint
        self.splrep = splrep
        self.splev  = splev
        self.splint = splint
        self.tck    = self.splrep(x,y,s=s)
    
    def __call__(self,x,der=0):
        """Evaluate function or its derivatives (der=1,2) at x."""
        return self.splev(x,self.tck,der=der)
        
    def integrate(self,x1,x2):
        """Integrate the function from x1 to x2."""
        return self.splint(x1,x2,self.tck)  
    
    
class Const:
    """ Simple class for constants. 
    
        Example:
            const=Const()
            E=const.kB*T
            const.list_constants() -> list of available constants
            const.search('Boltzmann')
    """   
    def __init__(self,lst):
        self.constants=[]
        const=lst.splitlines()
        for c in const:
            self.add_const(c)
    
    def add_const(self,st):
        """ Add constant (given as string) to list. """
        div=st.split(';')
        if len(div)<4: return
        for i in range(len(div)): div[i]=div[i].strip()
        div[2]=float(div[2])
        self.constants.append(div)
        exec('self.%s=%20.15g' %(div[1],div[2]) )
           
    def list_constants(self):
        """ List all the constants in this object. """
        for c in self.constants: print c
        
    def info(self,const):
        """ Print info for the given constant. """
        for c in self.constants:
            if c[1]==const.strip():
                print c
                return
        print 'Constant',const,'was not found.'
        
    def search(self,st):
        """ Search a part of string from the constant database. """
        for c in self.constants:
            hlp=[c[0],c[1],c[3]]
            for part in hlp:
                if part.find(st.strip())>=0:  
                    print c
                    exit
       
                    
SI_const="""
speed of light in a vacuum;     c;  2.99792458E+08;                 m/s 
elementary charge (of proton);  e;  1.60217733E-19;                 C
gravitational constant;         G;  6.67259E-11;                    m^3/(kg*s^2)
gravitational acceleration;     g;  9.80665;                        m/s^2
unified atomic mass constant;   u;  1.6605402E-27;                  kg
rest mass of electron;          me; 9.1093897E-31;                  kg
rest mass of proton;            mp; 1.6726231E-27;                  kg
rest mass of neutron;           mn; 1.6749286E-27;                  kg
Planck constant;                h;  6.6260755E-34;                  J*s
Planck constant/2pi;            hbar;   1.05457266E-34;             J*s
Hartree energy;                 Ha; 4.3597482E-18;                  J
Bohr radius;                    a_B;    5.29177249E-11;             m
Avogadro constant;              Na; 6.0221367E+23;                  1/mol
Faraday constant (Na*e);        F;  9.6485309E+04;                  C/mol
Boltzmann constant (R*Na);      kB; 1.380658E-23;                   J/K
Stefan-Boltzman constant;       sigma;  5.67051E-08;                J/(s*m^2*K^4)
Bohr magneton (e*hbar/(2*me));  mu_B;   9.2740154E-24;              J/T
nuclear magneton (e*hbar/(2*mp));   mu_N;   5.0507866E-27;          J/T
electric constant (1/(4*pi*epsilon_0)); Cc; 8987551787.37;          V*m/(A*s)
molar gas constant;             Rc; 8.314510;                       J/(mol*K)
permeability of a vacuum (4*pi*10E-7);  mu0;    1.25663706144E-7;   V*s/(A*m)   
permittivity of a vacuum;   epsilon_0;   8.85418781762E-12;          A*s/(V*m)
fine structure constant (e^2/(2*epsilon_0*h*c));    alpha;  7.29735308E-3; 1 """


AU_const="""
time unit;                      t_0;1.03274989966e-15;              s
speed of light in a vacuum;     c;  5850.79255585;                  L/T
elementary charge (of proton);  e;  1;                              Q
gravitational constant;         G;  7.97499949791e-37;              L^3/(M*T^2)
gravitational acceleration;     g;  1.97655923555e-19;              L/T^2
unified atomic mass constant;   u;  1;                              M
rest mass of electron;          me; 0.000548579895868;              M
rest mass of proton;            mp; 1.00727648749;                  M
rest mass of neutron;           mn; 11.0086648911;                  M
Planck constant;                h;  0.14716340209;                  E*T
Planck constant/2pi;            hbar;   0.0234217826822;            E*T
Hartree energy;                 Ha; 1;                              E
Bohr radius;                    a_B; 1;                             L
Avogadro constant;              Na; 6.0221367E+23;                  1/mol
Faraday constant (Na*e);        F;  9.6485309E+04;                  Q/mol
Boltzmann constant (R*Na);      kB; 3.1668297361e-06;               E/K
Stefan-Boltzman constant;       sigma;  3.76147527217e-26;          E/(T*L^2*K^4)
Bohr magneton (e*hbar/(2*me));  mu_B;   2.12719063966e-06;          E/Tesla
nuclear magneton (e*hbar/(2*mp));   mu_N;   1.15850422013e-09;      E/Tesla
permeability of a vacuum (4*pi*10E-7);  mu0;    3.67096691757e-08;   V*s/(A*m)   
permittivity of a vacuum;   epsilon_0;   0.0795774702667;           A*s/(V*m)
fine structure constant (e^2/(2*epsilon_0*h*c));    alpha;  7.29735308E-3; 1 """


SI=Const(SI_const)
AU=Const(AU_const)
         
        
if __name__=='__main__':
    pass
    #
    # testing class Const
    #
    #SI.list_constants()
    #SI.info('e')
    #SI.search('J')
    #SI.add_const('kk;test;3.4;XX')
    #SI.list_constants()
    from numpy import sqrt
    m=1/SI.a_B
    C=1/SI.e
    kg=1/SI.u
    J=1/( SI.me*SI.e**4/(8*SI.h**2*SI.epsilon_0**2)*2 )
    s=sqrt( kg*m**2/J )
    A=C/s
    
    print 'Basic quantities'
    print 'm',m
    print 'C',C
    print 'kg',kg
    print 'J',J
    print 's',s
    print 'A',A
    
    t_0     =1/s
    c       =SI.c*m/s
    G       =SI.G*m**3/(kg*s**2)
    g       =SI.g*m/s**2    
    me      =SI.me*kg
    mp      =SI.mp*kg
    mn      =SI.mn*kg
    h       =SI.h*J*s
    hbar    =SI.hbar*J*s
    Ha      =SI.Ha*J
    a_B     =SI.a_B*m
    F       =SI.F*C
    kB      =SI.kB*J
    sigma   =SI.sigma *J/(s*m**2)
    mu_B    =SI.mu_B*J
    mu_N    =SI.mu_N*J
    Cc      =SI.Cc *J*m/(C*A*s)
    mu0     =SI.mu0 *J*s/(A*m*C)
    epsilon_0=SI.epsilon_0 *A*s*C/(m*J)
    # C*V=J => V=J/C  
    
    print '\nConstants'
    print 't_0',t_0
    print 'c',c
    print 'G',G
    print 'g',g
    print 'me',me
    print 'mp',mp
    print 'mn',mn
    print 'h',h
    print 'hbar',hbar
    print 'Ha',Ha
    print 'a_B',a_B
    print 'F',F
    print 'kB',kB
    print 'sigma',sigma
    print 'mu_B',mu_B
    print 'mu_N',mu_N
    print 'Cc',Cc
    print 'mu0',mu0
    print 'epsilon_0',epsilon_0
    
    a=[1.0,2.0,3.22]
    print a2s(a)
       
       
    f=AnalyticFunction('x**2+a',a=2,b=5)
    f.info()
    print f(5,n=3,a=1)
    
              
 
        
        
        

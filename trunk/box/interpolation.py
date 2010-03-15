import numpy as nu
from box import mix
from scipy.linalg import norm
from scipy.optimize import fminbound
from scipy.optimize import brentq
from scipy.interpolate import splprep
from scipy.interpolate import splrep  
from scipy.interpolate import splev
from scipy.interpolate import splint
try:
    import pylab as pl
except:
    pass
vec=nu.array
linspace=nu.linspace

class MultipleSplineFunction:
    """ Class for interpolating many functions at the same time on the same grid. 
    
    cubic natural spline.
    """
    def __init__(self,x):
        """ Setup x-grid. 
        
        parameters:
        -----------
        x: the *linear* x-grid
        """
        self.x=x
        self.a=x[0]
        self.b=x[-1]
        self.n=len(x)
        self.xmin, self.xmax=x[0], x[-1]
        self.h=x[1]-x[0]
        self.labels=[]
        self.indices=[]
        self.y=[]        
        self.initialized=False
        self.m=0
        
    def add_function(self,y,label='',index=0):
        """ Add function on the same x-grid. 
        
        label: characterizes the function and 
        index: some additional number to characterize function (used for speed-up)
        """            
        assert not self.initialized
        self.y.append(y)
        self.labels.append(label)
        self.indices.append(index)
        self.m+=1
        
        
    def get_labels(self):
        return self.labels   
        
             
    def get_indices(self):
        return self.indices             
        
        
    def _initialize(self):
        """ Calculate second derivatives for all the functions. """            
        assert not self.initialized
        self.initialized=True

        self.d2=[]
        n=self.n
        for i in range(self.m):
            u=nu.zeros_like(self.x)
            d2=nu.zeros_like(self.x)
            y=self.y[i]
            x=self.x
            qn       = 0.0 
            un       = 0.0
    
            # Solve the second derivatives
            for i in range(1,n-1):
                sig  = (x[i]-x[i-1])/(x[i+1]-x[i-1])
                p    = sig*d2[i-1]+2
                d2[i]= (sig-1)/p
                u[i] = ( 6*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1] )/p
        
            d2[n-1]=(un-qn*u[n-2])/(qn*d2[n-2]+1)
            for i in range(n-1)[::-1]:
                d2[i]=d2[i]*d2[i+1]+u[i]
            self.d2.append(d2)
        self.d2=nu.array(self.d2)
        self.y=nu.array(self.y)            
    
                    
    def get_range(self):    
        return self.x[0],self.x[-1]
    
        
#    def _find_bin_fast(self,x):
#        """ For given x, return i such that x_i<=x<x_i+1 
#        Make it fast, but rely on linear scale.
#        """
#        lo=int( (x-self.xmin)/(self.xmax-self.xmin)*(self.n-1) )
#        return lo,lo+1
    
    
#    def _find_bin(self,x):
#        """ For given x, return i and i+1 such that x_i<=x<x_i+1 """
#        lo=0
#        hi=len(self.x)
#        while True:
#            if hi-lo>1:
#                k=(hi+lo)/2
#                if x<self.x[k]:
#                    hi=k
#                else:
#                    lo=k
#            else:
#                return lo, hi
#        raise AssertionError('Error in binding bin for %f' %x)
                  
  
    def __call__(self,x,der=None):
        """ Return interpolated values for all functions.
        
        parameters:
        x: x-point (If x is outside the interpolation range, return 0's)
        der: if None, return both function values and derivatives. If 0,1,2:
             return function values, first, or second derivatives.
        """
        if not self.initialized:
            self._initialize()
                  
        if x<self.x[0] or x>=self.x[-1]:
            if der==None:
                return nu.zeros((self.m,)),nu.zeros((self.m))
            else:
                return nu.zeros((self.m,))
               
        #lo, hi = self._find_bin(x)
        #xlo, xhi=self.x[lo], self.x[hi]
        hi = nu.searchsorted(self.x, x)
        lo = hi-1
        xlo, xhi = self.x[lo], self.x[hi]
        assert xlo<=x<=xhi
        
        h=self.h
        a, b=(xhi-x)/h, (x-xlo)/h
        ylo, yhi=self.y[:,lo], self.y[:,hi]
        dlo, dhi=self.d2[:,lo], self.d2[:,hi]
        
        y=a*ylo + b*yhi + ((a**3-a)*dlo+(b**3-b)*dhi)*(h**2)/6
        dy=(yhi-ylo)/h - (3*a**2-1)/6*h*dlo + (3*b**2-1)/6*h*dhi                
        
        if der==None:
            return y,dy
        elif der==0:
            return y
        elif der==1:
            return dy
        elif der==2:
            return a*self.d2[:,lo] + b*self.d2[:,hi]


class FastSplineFunction:
    def __init__(self,x,y,grid,**kwargs):
        """ Initialize second derivatives.
        
        Parameters:
        -----------
        x: x-grid 
        y: y-values on given x-grid
        grid: type of grid
            'linear':   x[i]=xmin + i/(N-1)*(xmax-xmin) (i~(N-1)*(x-xmin)/(xmax-xmin))
            'power':  x[i]=xmin + (i/(N-1))**p*(xmax-xmin)
        **kwargs: parameters for grid (can vary) (xmin,xmax,p,...) 
            (N comes from length of x)
        """
                
                
        self.x=x
        self.y=y
        self.a=x[0]
        self.b=x[-1]
        d2=nu.zeros_like(y)
        u=nu.zeros_like(y)
        n=len(x)
        qn       = 0.0 
        un       = 0.0
  
        # Solve the second derivatives
        for i in range(1,n-1):
            sig  = (x[i]-x[i-1])/(x[i+1]-x[i-1])
            p    = sig*d2[i-1]+2
            d2[i]= (sig-1)/p
            u[i] = ( 6*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1] )/p
    
        d2[n-1]=(un-qn*u[n-2])/(qn*d2[n-2]+1)
        for i in range(n-1)[::-1]:
            d2[i]=d2[i]*d2[i+1]+u[i]
        self.n=n
        self.d2=d2
        self.grid=grid
        self.xmin, self.xmax=x[0], x[-1]
        if grid=='linear':
            #self.xmin,self.xmax=kwargs['xmin'],kwargs['xmax']
            self.find_lower_index=self.find_lower_index_linear
        if grid=='power':
            self.p=kwargs['p']
            self.find_lower_index=self.find_lower_index_power
            
    def get_range(self):    
        return self.x[0],self.x[-1]
        
    def find_lower_index_linear(self,x):
        """ For given x, return i such that x_i<=x<x_i+1 """
        return int( (x-self.xmin)/(self.xmax-self.xmin)*(self.n-1) )
        
    def find_lower_index_power(self,x):
        """ For given x, return i such that x_i<=x<x_i+1 """
        return int( ((x-self.xmin)/(self.xmax-self.xmin))**(1.0/self.p)*(self.n-1) )
  
    def __call__(self,x,der=0):
        """! If x is outside the interpolation range, return 0."""
                  
        if x<self.x[0] or x>=self.x[-1]:
            return 0.0
        
        lo=self.find_lower_index(x)
        hi=lo+1
        
        h=self.x[hi]-self.x[lo]
        assert self.x[lo]<=x<self.x[hi]
        a, b=(self.x[hi]-x)/h, (x-self.x[lo])/h
  
        if der==0:
            return a*self.y[lo] + b*self.y[hi] + ((a**3-a)*self.d2[lo]+(b**3-b)*self.d2[hi])*(h**2)/6
        elif der==1:
            return (self.y[hi]-self.y[lo])/h - (3*a**2-1)/6*h*self.d2[lo] + (3*b**2-1)/6*h*self.d2[hi]
        elif der==2:
            return a*self.d2[lo] + b*self.d2[hi]
  


class SplineFunction:
    def __init__(self,x,y,k=3,s=0,name=None):
        """ Simple B-spline function; order k is cubic by default. 
        
        Parameters:
        -----------
        x:  x-grid
        y:  y-values for given x-grid points
        k:  order of spline (cubic by default)
        s:  smoothness parameters (means exact reproduction of x,y values)
        name: name of the function
        
        """        
        if s==-1:
            s=len(x)-nu.sqrt(2.0*len(x))
        self.tck=splrep(x,y,s=s,k=k)
        self.x=x
        self.y=y
        self.a=x[0]
        self.b=x[-1]
        self.M=len(y)
        self.name=name
        
    def __call__(self,x,der=0):
        """ Return der'th derivative of f(x) 
        
        Return zero if x beyond the original grid range.
        """
        if x<self.x[0] or self.x[-1]<x:
            return 0.0
        else:
            return splev(x,self.tck,der=der)
        
    def get_name(self):
        """ Return the name of the function. """
        return self.name        
    
    def get_range(self):
        return (self.x[0],self.x[-1])
    
    def limits(self):
        return self.get_range()
    
    def solve(self,y,a=None,b=None):
        """ Solve x for f(x)=y, where x in [a,b]. """
        if a==None: a=self.x[0]
        if b==None: b=self.x[-1]
        assert a<b
        return brentq(lambda x:self(x)-y,a=a,b=b)
        
    def integrate(self,a,b):
        """
        Integrate given function within [a,b]
        """
        return splint(a,b,self.tck)
    
    def plot(self,return_pylab=False,der=0,filename=None):
        """
        Plot the function with matplolib
        """
        p=pl
        p.clf()
        from numpy import linspace
        X=linspace(self.a,self.b,self.M*10)
        Y=[self(x,der=der) for x in X]
        p.plot(X,Y)
        if der==0:
            p.scatter(self.x,self.y)
        f1=SplineFunction(self.x,self.y,k=1)
        p.plot(X,[f1(x,der=der) for x in X])
        if return_pylab:
            return p
        elif filename!=None:
            p.savefig(filename)
        else:
            p.show()
        p.close()
    
    def max_deviation_from_linear(self):
        """
        For given spline (default cubic), return maximum difference
        wrt linear (k=0) interpolation.
        """
        from numpy import array,linspace,abs 
        min_dx=min( array(self.x[1:])-array(self.x[0:-1]) )
        M=10*(self.b-self.a)/min_dx
        X=linspace(self.a,self.b,M)
        f1=SplineFunction(self.x,self.y,k=1)     
        return max( array([abs(f1(x)-self(x)) for x in X]) )
                
    def smoothness(self):
        """
        Return a measure for the ''smoothness'' of interpolation.
        
        It is measured by max_deviation_from_linear/average|dy|.
        Smooth interpolations should have <<1.
        If ~1, then interpolation deviates as much as is the variations
        in y and interpolation is not smooth.        
        """
        from numpy import abs,average,array
        avg=average( abs(array(self.y[1:])-array(self.y[0:-1])) )        
        return self.max_deviation_from_linear()/avg


class Function:
    
    def __init__(self,mode,*args,**kwargs):
        if mode is 'spline':
            self.f=SplineFunction(*args,**kwargs)
        elif mode is 'string':
            raise NotImplementedError('todo')
            self.args=args
            self.kwargs=kwargs
        elif mode is 'fastspline':
            self.f=FastSplineFunction(*args,**kwargs)
        else:
            raise NotImplementedError('todo')
        
    def __call__(self,x,der=0):
        return self.f(x,der)
        
    def plot(self,der=0,a=None,b=None,npoints=1000,filename=None,return_pylab=False):
        """ Plot the function with matplolib. """
        import pylab as pl
        from numpy import linspace
        a0,b0=self.f.get_range()
        lower=[a0,a][a!=None]
        upper=[b0,b][b!=None]
        X=linspace(lower,upper,npoints)
        Y=[self(x,der=der) for x in X]
        pl.plot(X,Y)
        if return_pylab:
            return pl
        elif filename!=None:
            pl.savefig(filename)
        else:
            pl.show()
        pl.close() 

#class Function(SplineFunction,FastSplineFunction):
    
    #def __init__(self,mode,*args,**kwargs):
        #if mode is 'spline':
            #SplineFunction.__init__(self,*args,**kwargs)
            #self.__call__=SplineFunction.__call__
        #elif mode is 'string':
            #raise NotImplementedError('todo')
            #self.args=args
            #self.kwargs=kwargs
        #elif mode is 'fastspline':
            #FastSplineFunction.__init__(self,*args,**kwargs)
            #self.__call__=FastSplineFunction.__call__
        #else:
            #raise NotImplementedError('todo')
        
    #def plot(self,der=0,a=None,b=None,npoints=1000,filename=None,return_pylab=False):
        #""" Plot the function with matplolib. """
        #import pylab as pl
        #from numpy import linspace
        #lower=[self.a,a][a!=None]
        #upper=[self.b,b][b!=None]
        #X=linspace(lower,upper,npoints)
        #Y=[self(x,der=der) for x in X]
        #pl.plot(X,Y)
        #if return_pylab:
            #return pl
        #elif filename!=None:
            #pl.savefig(filename)
        #else:
            #pl.show()
        #pl.close()   
        

def grid(n):
    return nu.linspace(0,1,n)



class VectorSplineFunction:
    def __init__(self,r,k=3,s=None,u=None,force=None):
        """ 
        Interpolate vector r's path in N-dimensional space (from M points {r}).
        
        r=array([[vec1], [vec2], [vec3], ..., [vecM]])
        (r[:,i] is the trajectory of component i.)
        
        k = order of the spline (cubic (3) by default to get curvature)
        s = smoothing parameter ala scipy (s=0 -> points r given exactly)
        u = node points (for internal u-parameter)
        force = force field at the nodes of r (forces evaluated at each M point)
        
        Usage:
            fr=VectorSplineFunction(r) (now r=fr(t) with t [0,1])
            tan=fr.normalized_tangent(0.2)
            ...            
            
        
        The trajectory is parametrized with variable 't' ([0,1]) where 
        t=l/L (L is the total length of the trajectory and l length so far).
        Internally for technical reasons the path is internally parametrized 
        by 'u' ([0,1]), where the path is taken to be _linear_ between 
        the points. Thus, for k=1, u is equal to t.   
        
        r=r(t) and r=r(u); r(t=0)=r(u=0) and r(t=1)=r(u=1)
        The length along the path l(u) and l(t). t has the property
        dl(t)/dr=L/1=L = constant (l(t=0)=0 and l(t=1)=L).
        
        dl   dl(t) dt   dt         l'(u)
        -- = -----*--=L*--  --> dt=-----du --> t(u)=l(u)/L (with correct bound.cond.)
        du    dt   du   du           L
        
                 /b       /b
        l(a,b)= |   ds = |    |dr(u)/du|du          (for a,b=u)
                /u=a     /u=a
                
                 /b       /b
        l(a,b)= |   ds = |    |dr(u)/du*du/dt|dt    (for a,b=t)
                /t=a     /t=a
                
        where |...| is the Eucledian norm.
        
        
        Energy slope E'(t)=F(t)*r'(t)

        """
        self.k=k
        self.s=s
        self.r=r.copy()
        self.u=u
        self.N=len(r[0,:])
        self.M=len(r[:,0])
        self.M2=self.M*10 #finer grid  (used for plotting etc)
        self.M3=self.M*50 #finest grid (used in integration etc) 
        self.w=[1]*self.M
        self.w[0]=1E6
        self.w[-1]=1E6
        
        if self.s==None:
            self.s=self.M-nu.sqrt(2.0*self.M)
        if self.u!=None:
            self.u=u.copy()
        
        self.tck,self.u=splprep(self.r.transpose(),u=self.u,k=self.k,s=self.s,w=self.w)            
        
        u,l=self.line_integral(a=0,b=1,parameter='u',full_out=True)
        self.length=l[-1]
        t=l/self.length
        self.u_of_t=SplineFunction(t,u,k=3,s=0)
        self.t_of_u=SplineFunction(u,t,k=3,s=0)
        self.t=[self.t_of_u(u) for u in self.u]
         
        if self.k!=1: 
            self.linear=VectorSplineFunction(self.r,k=1,s=0,u=self.u)
        if force!=None:
            self.force=VectorSplineFunction(r=force,k=self.k,s=self.s,u=self.u)
       
        
    def get_length(self):
        return self.length
                       
                               

                               
    def plot(self,out='screen'):
        selected=range(self.N) #from maximum deviation?
        n=len(selected)
        ny=n+2
        nx=1
        
        rhomo=self.homogenize()
        thomo=grid(self.M)
        tsmooth=grid(self.M2)
        rsmooth=vec([self(t=t) for t in tsmooth])
        
        
        fig=1
        if self.N==2:
            pl.figure(fig)
            pl.scatter(self.r[:,0],self.r[:,1],color='r',label='trajectory')
            pl.scatter(rhomo[:,0],rhomo[:,1],s=[50]*self.M,color='y')
            pl.plot(rsmooth[:,0],rsmooth[:,1])            
            fig+=1
        
        # plot u(t)
        pl.figure(fig)
        pl.subplot(ny,nx,1)
        usmooth=[self.u_of_t(x) for x in tsmooth]
        pl.plot(tsmooth,usmooth)
        pl.plot(tsmooth,tsmooth)
        pl.title('u(t)')
        pl.axis(xmin=0,xmax=1.0,ymin=0,ymax=1)
        
        pl.subplot(ny,nx,2)
        pl.plot(tsmooth,self.curvature(t=tsmooth),label='curvature')
        
        sub=3
        for i in selected:
            pl.subplot(ny,nx,sub)
            pl.plot(tsmooth,rsmooth[:,i])
            pl.scatter(thomo,rhomo[:,i],s=[100]*self.M,label='homog',color='y',marker='d')
            pl.scatter( self.t,self.r[:,i],color='r',label='%i' %i)
            pl.axis(xmin=0,xmax=1.0) 
            sub+=1
            #p.legend()
            
        if out=='screen':
            pl.show()
        else:
            pl.savefig(out)
            pl.close()
        
            
                    
    def __call__(self,t=None,u=None,der=0):
        """ 
        Return the interpolated r(t) or r(u). If you evaluate t-derivatives
        then chain rule is used: dr(t)/dt=dr(u)/du*(du/dt). This because
        spline routine is always called using internal u-parameter. Take special
        care to return initial or end-points if parameter==0 or 1.
        """
        assert t!=None and u==None or t==None and u!=None
        if t!=None:
            u=self.u_of_t(t)
        elif u!=None:
            pass
        
        if der==1 and t!=None:
            return vec( splev(u,self.tck,der=1) )*self.u_of_t(t,der=1) 
        elif der==2 and t!=None:
            return vec(splev(u,self.tck,der=2))*self.u_of_t(t,der=1)**2 + \
                   vec(splev(u,self.tck,der=1))*self.u_of_t(t,der=2)
        else:
            if abs(u)<1E-4 and der==0:
                return self.r[0,:]
            elif abs(u-1)<1E-4 and der==0:
                return self.r[-1,:]
            else:
                return vec(splev(u,self.tck,der=der))
        
        
    def line_integral(self,a=None,b=None,parameter=None,full_out=False):
        """ 
        Return line integral in N-dimensional space along splined function.
        Paramater is 'u' or 't' [a,b]. If b<a, the integral is negative.
        if full_out==True, return the function where the upper limit goes
        from a to b.
        """
        div=self.M3 #use finest grid in intergration
        sign=1
        if b<a:
            a,b=(b,a)
            sign=-1
        
        dx=(b-a)/(div-1.0)
        grid=nu.linspace(a,b,div) 
        
        if parameter=='u':
            ders=vec([norm(self(u=x,der=1)) for x in grid])
        elif parameter=='t':
            ders=vec([norm(self(t=x,der=1)) for x in grid])
            
        l=[0.0] # upper limit=a
        for i in range(len(grid)-1):
            l.append(l[-1]+ders[i]+ders[i+1])
        l=vec(l)/2*dx*sign
        
        if full_out:
            return vec(grid),vec(l)
        else:
            return l[-1]
    
    def tangents(self,t=None):
        """ Return tangents at t (list). Nodes by default. """
        if t==None:
            t=self.t
        return [self.tangent(tx) for tx in t]
            
    def tangent(self,t):
        """ Return tangent r'(t) """
        return self(t=t,der=1)
            
    def normalized_tangent(self,t):
        """ Return r'(t)/|r'(t)| """
        tangent=self.tangent(t)
        return tangent/norm(tangent)
    
    def closest_parameter(self,r):
        """ 
        Return the parameter t that minimizes 
        |r(t)-r| along the splined curve. 
        """        
        t=grid(self.M2)
        d=[norm(self(t=tx)-r) for tx in t]
        i=nu.argmin(d)
        tx=fminbound(lambda tx:norm(self(t=tx)-r),t[max(0,i-2)],t[min(self.M2-1,i+2)])
        return tx
        
    def closest_parameters(self,r=None):
        """ Return the parameters t (in [0,1]) minimizing |f(t)-r| for all r. """        
        if r==None:
            r=self.r
        t=[self.closest_parameter(rx) for rx in r]
        return t
        
    def homogenize(self,div=None):
        """ Return set of r that have equal distances in between. """
        if div==None:
            div=self.M
        t=grid(div)
        return vec([self(t=tx) for tx in t])
    
    def deviation_from_linear(self,t):
        """ Return |r(t)-rl(t')| where rl is linear and t' minimizing the norm."""
        rt=self(t=t)
        rl=self.linear(t=self.linear.closest_parameter(rt))
        return norm(r-rl)
        
    def max_deviation_from_linear(self):
        """ Return the maximum deviation from linear interpolation. """
        deviation=self.deviation_from_linear
        t=grid(self.M2)
        dev=[deviation(tx) for tx in t]
        i=nu.argmax(dev)
        tmax=fminbound(lambda tx:-deviation(tx),t[max(0,i-2)],t[min(self.M2-1,i+2)])
        return deviation(tmax)
        
    def deviations_at_nodes(self):
        """ Return |r(t_i)-r_i| for all nodes. """
        raise NotImplementedError()
    
    def curvature(self,t):
        """ Return the curvature norm |r''(t)| at t (can be a list)."""
        #assert self.k>2 # at least cubic required
        if type(t)==type(1.0) or type(t)==type(1):
            return norm(self(t=t,der=2))
        else:
            return [norm(self(t=tx,der=2)) for tx in t]
        

class BilinearInterpolation:
    """ 
    Perform bilinear interpolation for 2D (xy)-data
    
    For bilinear interpolation, see e.g. Wikipedia.
    """
    def __init__(self,a,xgrid=None,ygrid=None):
        self.a=a
        self.nx=a.shape[0]
        self.ny=a.shape[1]
        self.xgrid=xgrid
        self.ygrid=ygrid
        if self.xgrid==None:
            self.xgrid=nu.arange(self.nx)
        if self.ygrid==None:
            self.ygrid=nu.arange(self.ny)
        self.dx=self.xgrid[1]-self.xgrid[0]
        self.dy=self.ygrid[1]-self.ygrid[0]
        
    def __call__(self,x,y):
        """ Return the interpolated value at x,y. """
        assert self.xgrid[0]<=x<=self.xgrid[-1]
        assert self.ygrid[0]<=y<=self.ygrid[-1]                
        
        i,j=(nu.floor(x/self.dx),nu.floor(y/self.dx))
        i,j=(min(i,self.nx-2),min(j,self.ny-2))
        dx,dy=((x%self.dx)/self.dx,(y%self.dy)/self.dy)
        raise NotImplementedError('Check that dx=0...1 and not something else.')
        
        a=self.a
        return a[i,j] *(1-dx)*(1-dy)+a[i+1,j]*dx*(1-dy)+\
               a[i,j+1]*(1-dx)*dy   +a[i+1,j+1]*dx*dy
        
    def plot(self,nx=100,ny=100,out='screen'):
        """ Make 2D contour-color plot of the data. """
        b=nu.empty((nx,ny))
        X=nu.linspace(self.xgrid[0],self.xgrid[1],nx)
        Y=nu.linspace(self.ygrid[0],self.ygrid[1],ny)
        for i,x in enumerate(X):
            for j,y in enumerate(Y):
                b[i,j]=self(x,y)
        pl.contourf(b,100)
        pl.hot()
        pl.colorbar()
        if out=='screen':
            pl.show()
        else:
            pl.savefig('linter.png')
    
class TrilinearInterpolation:
    """ 
    Perform trilinear interpolation for 3D (xyz)-data
    
    For trilinear interpolation, see e.g. Wikipedia.
    """
    def __init__(self,a,grids=None):
        self.a=a
        self.n=a.shape
        self.grids=grids
        if self.grids==None:
            self.grids=vec([nu.arange(self.n[i]) for i in range(3)])*1.0
        self.dg=vec([self.grids[i][1]-self.grids[i][0] for i in range(3)])*1.0
                
    def __call__(self,r):
        """ Return the interpolated value at r. """
        for i in range(3):
            assert self.grids[i][0]<=r[i]<=self.grids[i][-1]       

        ind0=nu.floor((r/self.dg))
        dx0,dy0,dz0=[self.grids[i][1]-self.grids[i][0] for i in range(3)]
        # i,j,k grid point can never be a point at 'larger' side
        i,j,k=[int(min(ind0[i],self.n[i]-2)) for i in range(3)]
        dx,dy,dz=[r[p]-self.grids[p][ind] for p,ind in zip(range(3),(i,j,k))]
        dx,dy,dz=(dx/dx0, dy/dy0, dz/dz0)
                            
        a=self.a
        i1=a[i,j,k]    *(1-dz)+a[i,j,k+1]*dz
        i2=a[i,j+1,k]  *(1-dz)+a[i,j+1,k+1]*dz
        j1=a[i+1,j,k]  *(1-dz)+a[i+1,j,k+1]*dz
        j2=a[i+1,j+1,k]*(1-dz)+a[i+1,j+1,k+1]*dz
        w1=i1*(1-dy)+i2*dy
        w2=j1*(1-dy)+j2*dy
        return w1*(1-dx)+w2*dx     

        
    def write_vtk(self,name='trilinear',gpts=None):
        """ 
        Write data into vtk file with given number of grid points. 
        
        Parameters:
        gpts - number of grid points to all directions. Default
               is the number of DFT grid points.
        """
        if gpts==None:
            gpts=self.n
        R=[nu.linspace(0,self.grids[i][-1],gpts[i]) for i in range(3)]
        nx,ny,nz=(gpts[0],gpts[1],gpts[2])
        
        
        of=open('%s.vtk' %name,'w')
        print>>of,"# vtk DataFile Version 2.0"
        print>>of,"Trilinear interpolation"
        print>>of,"ASCII"
        print>>of,"DATASET RECTILINEAR_GRID"
        print>>of,"DIMENSIONS %i %i %i" %(nx,ny,nz)
        print>>of,"X_COORDINATES %i double" %nx
        print>>of,mix.a2s(R[0])
        print>>of,"Y_COORDINATES %i double" %ny
        print>>of,mix.a2s(R[1])
        print>>of,"Z_COORDINATES %i double" %nz
        print>>of,mix.a2s(R[2])
        print>>of,"POINT_DATA %i" %(nx*ny*nz)
        print>>of,"SCALARS data double"
        print>>of,"LOOKUP_TABLE default"
        
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    r=vec([R[0][i],R[1][j],R[2][k]])
                    print>>of, self(r)
        of.close()
        
     
        
if __name__=='__main__':
    #sp=SplineFunction([0,1,2,3],[0,1,2,10])
    #print sp.solve(4.56)
    #sp.plot()
    #assert False
    if False:
        r=vec([vec([1,0]),vec([0.5,1.0]),vec([0,2]),vec([0,2.1]),vec([0.2,2.2]),vec([0.2,2.3]),vec([0.21,2.2]),vec([0.2,2.3]),vec([-0.3,2.5])])
        #r=vec([[0,0],[1,1],[2,10],[3,8]])
        f=VectorSplineFunction(r,k=3,s=0)
        #for t in nu.linspace(0,1,10):
            #print t,f(t,der=0) #,f(t,der=1),f(t,der=-1)
            #print t, f.line_integral(b=t) #,nu.sqrt(5.0)
        f.max_deviation_from_linear()
        f.plot()
        #print 'length',f.get_length()
        
    if False:
        a=vec([[10,1],[2,5]])
        inte=BilinearInterpolation(a)
        inte.plot()
        
    a=nu.empty((2,2,2))
    a[0,0,0]=0.1
    a[0,0,1]=1
    a[0,1,0]=2
    a[0,1,1]=3
    a[1,0,0]=4
    a[1,0,1]=4
    a[1,1,0]=3
    a[1,1,1]=4
    inte=TrilinearInterpolation(a)
    print inte(vec([0.999,0.999,0.999]))
    print inte(vec([1.0,1.0,1.0]))
    inte.write_vtk(gpts=[20,20,20])

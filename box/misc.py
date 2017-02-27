import numpy as np
from scipy.linalg import solve
from box.dummymixer import DummyMixer

class AndersonMixer(DummyMixer):
    """
    produces the Anderson mixing of input vectors in iterative process.
    See, e.g. V. EYERT in J. Comp. Phys. 124, 271 (1996) (the notation
    is almost directly from there)
    x is the input, and y the output vector of iteration. F is the
    history of residuals. Mmax is the maximum number of last iterations
    taken into account.


    """
    def __init__(self, mixing_constant=0.2, memory=3, convergence=1e-3, chop=None):
        DummyMixer.__init__(self, mixing_constant, convergence)
        self.name = 'Anderson'
        self.memory = memory
        self.chop=chop

    def set(self,**kwargs):
        """ Reset parameters given in initialization. """
        if 'memory' in kwargs:
            raise KeyError('Not allowed to change memory on-fly.')
        self.__dict__.update(**kwargs)

    def __call__(self,xi,yi):
        if self.it==0:
            self.F=np.zeros((self.memory+1,len(xi)))
            self.x=np.zeros((self.memory+1,len(xi)))
        self.it+=1
        M=min(self.it-1,self.memory) #use M _previous_ iterations (it-1,it-2,...)
        self.F[0,:]=yi-xi
        self.x[0,:]=xi.copy()
        if M<=1:
            simple=True
        elif M>1:
            A=np.zeros((M,M))
            b=np.zeros((M,))
            # prepare for solving solve A*z=b (eq.(4.3) in Eyert)
            for i in range(M):
                b[i]=np.dot(self.F[0,:]-self.F[i+1,:],self.F[0,:])
                for j in range(M):
                    A[i,j]=np.dot(self.F[0,:]-self.F[i+1,:],self.F[0,:]-self.F[j+1,:])
            try:
                # c is the optimum linear combination
                c=solve(A,b)
                simple=False
            except:
                # if A singular, use simple mixing
                simple=True

        if simple:
            xb=(1-self.beta)*xi+self.beta*yi
        else:
            # xb is the best (calculated) estimate; go that direction
            xb=xi.copy()
            Fb=self.F[0,:].copy()
            for j in range(M):
                xb+=c[j]*(self.x[j+1,:]-xi)
                Fb+=c[j]*(self.F[j+1,:]-self.F[0,:])
            xb=xb+self.beta*Fb


        # The input must not change more than chop for all
        maxdev=max(abs(xb-xi))
        if self.chop!=None and maxdev>self.chop:
            xb = xi + self.chop/maxdev*(xb-xi)

        # shift history
        self.F[1:self.memory+1,:] = self.F[0:self.memory,:]
        self.x[1:self.memory+1,:] = self.x[0:self.memory,:]

        fmax=abs(self.F[0,:]).max()
        self.fmax.append(fmax)
        if fmax<self.convergence:
            done=True
        else:
            done=False
        return done, xb


    def out_of_iterations(self,out):
        """ Print out-of-iterations info. """
        print('Anderson mixer out of iterations!', file=out)
        print('History of maximum residuals:', file=out)
        for it,fmax in enumerate(self.fmax):
            print(it,fmax, file=out)

        print("Recent history: (number=i'th last iteration)", file=out)
        for i in range(1,self.memory+1):
            print(i-1, self.x[i,:5], file=out)

if __name__=='__main__':
    def testf(x):
        return -np.array([3.0,2.0,1.5,1.0,0.5])*x-0.01*x**3

    mix=AndersonMixer(0.2,3,1E-40)
    xi=np.array([1.,1.,1.,1.,1.])
    for i in range(50):
        out=testf(xi)
        #print xi,out
        done,xi=mix(xi,out)
        mix.echo()

        if done:
            break

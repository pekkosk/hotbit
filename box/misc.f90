class AndersonMixer:
    """
    produces the Anderson mixing of input vectors in iterative process.
    See, e.g. V. EYERT in J. Comp. Phys. 124, 271 (1996) (the notation
    is almost directly from there)
    x is the input, and y the output vector of iteration. F is the 
    history of residuals. Mmax is the maximum number of last iterations
    taken into account.
    """
    def __init__(self,beta,mmax,limit):
        self.mmax=mmax
        self.beta=beta
        self.limit=limit
        self.it=0
        self.fmax=[]
        
    def __call__(self,xi,yi):
        if self.it==0:
            self.F=nu.zeros((self.mmax+1,len(xi)))
            self.x=nu.zeros((self.mmax+1,len(xi)))
        self.it+=1
        M=min(self.it-1,self.mmax) #use M _previous_ iterations (it-1,it-2,...)
        self.F[0,:]=yi-xi
        self.x[0,:]=xi.copy()
        if M<=1:
            simple=True
        elif M>1:
            A=nu.zeros((M,M))
            b=nu.zeros((M,))
            # prepare for solving solve A*z=b (eq.(4.3) in Eyert)
            for i in range(M):
                b[i]=nu.dot(self.F[0,:]-self.F[i+1,:],self.F[0,:])
                for j in range(M):
                    A[i,j]=nu.dot(self.F[0,:]-self.F[i+1,:],self.F[0,:]-self.F[j+1,:])
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
            #print '...',xb,Fb,c
            #print M,self.it
            #print 'A',A,b
            #print '0',self.F[0,:]
            #print '1',self.F[1,:]
            #print '2',self.F[2,:]
            
        # shift history
        self.F[1:self.mmax+1,:] = self.F[0:self.mmax,:]
        self.x[1:self.mmax+1,:] = self.x[0:self.mmax,:]
        
        fmax=abs(self.F[0,:]).max()
        self.fmax.append(fmax)
        if fmax<self.limit:
            done=True
        else:
            done=False
        return done, xb

    def echo(self):
        """ Say something about the progress of iteration. """
        print 'Anderson:',self.fmax[-1]
        
    def final_echo(self):
        print self.it

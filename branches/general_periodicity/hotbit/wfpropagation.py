import numpy as nu
from box import mix
from hotbit.auxil import wf_dot_product, matrix_has_property, same_matrix
dot=nu.dot
hbar=0.02342178268

def matrix_square_root(A):
    """ Return the square root of real symmetric matrix A.
    
    Matrix multiplications to get inv, sqrt or inverse sqrt
    A_inv=U*sqrt(D)*U', where D are diagonal eigenvalue matrix and
    U'=transpose(U)
    """
    assert matrix_has_property(A,'symmetric')
    D, U=nu.linalg.eigh(A)
    assert all(D>=0)
    assert matrix_has_property(U,'orthogonal')
    d=nu.sqrt(D)
    return dot( dot(U,nu.diag(d)),U.transpose().conjugate() )


def density_matrix(wf,occ):
    n=len(occ)
    raise NotImplementedError('complexity?')
    rho=nu.zeros((n,n))
    for i in range(n):
        for j in range(n):
            rho[i,j]=sum( occ*wf[i,:]*wf[j,:].conjugate() )
    return rho
                    
                    
def mulliken(rho,calc,S):
    """ Return excess Mulliken populations. """
    q=[]
    N=len(calc.el)
    rhoS=dot(rho,S)
    for i in range(N):
        orbitals=calc.el.orbitals(i,indices=True)
        q.append( sum(rhoS.diagonal()[orbitals]) )           
    return nu.array(q)-calc.el.get_valences()


class WFPropagation:
    def __init__(self,calc):
        self.calc=calc
        self.st=calc.st
        self.ia=calc.ia
        self.es=calc.es
        self.env=calc.env
        
        # this is just in case ions are static...
        self.S=self.ia.S
        self.H0=self.ia.H0
        self.Si=nu.linalg.inv(self.S)
        self.Ss=matrix_square_root(self.S)
        self.Sis=nu.linalg.inv(self.Ss)
        self.n=self.S.shape[0]
        
        
    def propagate(self,dt):
        """                             -1/2             -1/2      -1/2           -1/2
                      -1/2    1-i/4*[ S(1)  * H(1) * S(1)     + S(0)  * H(0) * S(0)   ]*dt       1/2        
        wf(t+dt) = S(1)  *   --------------------------------------------------------------- S(0)    wf(t)
                                        -1/2             -1/2      -1/2           -1/2
                              1+i/4*[ S(1)  * H(1) * S(1)     + S(0)  * H(0) * S(0)   ]*dt               
        """ 
        raise RuntimeError('WF propagation still work in progress -- to be developed.')
        wf=self.st.wf
        if not isinstance(wf[0,0],complex):
            wf=wf*(1.0+0.0j) #make wave functions of complex type
        
        S=self.S
        H0=self.H0
        Si=self.Si # !!! calculate inverse, square root, and inverse square root
        Ss=self.Ss # !!! again if ions are moving
        Sis=self.Sis
        one = nu.identity(self.n)
        
        self.es.construct_tables()
        rho=density_matrix(wf,self.calc.st.f)
        dq=mulliken(rho,self.calc,S)
        
        H1=self.es.construct_H1( dq )
        H = H0 + H1*S
                
        SHS1=dot(Sis,dot(H,Sis))
        SHS0=SHS1    
        num = one - 0.25j*(SHS1+SHS0)*dt/hbar 
        den = one + 0.25j*(SHS1+SHS0)*dt/hbar 
        prod = dot( num,nu.linalg.inv(den) )
        
        operator = dot( Sis,dot(prod,Ss) )   
                       
        # check that matrices are OK; to be removed when thigs work (for efficiency)                       
        tol=1E-10
        assert same_matrix( nu.dot(Si,S),nu.identity(self.n),tol=tol )
        assert same_matrix( nu.dot(Ss,Ss),S,tol=tol )
        assert same_matrix( nu.dot(Sis,Sis),Si,tol=tol )
        assert same_matrix( nu.dot(Sis,Ss),nu.identity(self.n),tol=tol )
        
        # perform wave function propagation
        for i in range(self.n):
            wf[:,i]=dot(operator,wf[:,i])
                
        rho=density_matrix(wf,self.calc.st.f)
        dq=mulliken(rho,self.calc,S)        
        self.st.wf=wf   
        self.st.dq=dq  
        # propagation would need to be done self-consistenly (wrt. dq)
        
        # !!! make update to work via states.update
        #self.st.update(self.st.e,wf)
        self.env.propagate_time(dt)
        self.st.large_update()        
            
            
                    
                    
            
        
        
        

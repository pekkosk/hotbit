import numpy as np


def k_to_kappa_points(k,atoms):
    """ Transform normal k-points into kappa-points.
    
    parameters:
    ===========
    kpts:     list of k-vectors
    atoms:    hotbit.Atoms -object
    """
    if atoms.container.type!='Bravais':
        raise AssertionError('For other than Bravais lattices only kappa-points are allowed.')
    L = atoms.get_cell()
    kappas = []
    for kpt in k:
        kappas.append( np.dot(L,kpt) )
    return np.array(kappas)


def separate_symbols(string):
    """ From given string, separate element symbols AuC=>Au,C"""
    if string[1].islower():
        return string[:2],string[2:]
    else:
        return string[0],string[1:]
    
    
def make_hermitian(H):
    """ Make matrix H hermitian using fast numpy routines. """
    H2=H+H.transpose().conjugate()
    return H2-H.diagonal()*np.identity(H.shape[0])
        
        
def make_derivative_antihermitian(dH):
    """ Make vector matrix H antihermitian using fast numpy routines. """
    dH2=np.zeros_like(dH)
    for a in range(3):
        dH2[:,:,a]=dH[:,:,a]-dH[:,:,a].transpose().conjugate()
    return dH2
    
def wf_dot_product(wf1,wf2,S):
    """ Return the inner product of two wave functions wrt. metric S. """
    return np.dot( wf1.conjugate(),np.dot(S,wf2) )
    
    
def same_matrix(M1,M2,tol=1E-13): 
    return np.all( abs(M1.flatten()-M2.flatten())<tol )    
    
    
def matrix_has_property(A,property,tol=1E-12):
    """ Check whether matrix A has given property.
        
    Parameters:
    -----------
    A: General matrix
    property: property can be one of:
            'square' 
            'hermitian' 
            'symmetric'
            'orthogonal'
            'unitary'
    """
    property=property.lower()
    if property=='square':
        if len(A.shape)>2: return False
        return A.shape[0]==A.shape[1]
    elif property=='hermitian':
        return same_matrix(A, A.conjugate().transpose(), tol )
    elif property=='symmetric':                            
        return same_matrix(A, A.transpose(), tol )
    elif property=='orthogonal':                            
        return same_matrix( np.dot(A,A.transpose()), np.identity(A.shape[0]), tol )
    elif property=='unitary':
        #return same_matrix( np.dot(A,np.linalg.inv(A)), np.identity(A.shape[0]), tol )
        return same_matrix( np.dot(A,A.conjugate().transpose()), np.identity(A.shape[0]), tol )
    else:
        raise NotImplementedError('Property %s is not defined yet.' %property)
    
    
if __name__=='__main__':
    A=np.array([[np.cos(1.3),np.sin(1.3)],\
                [-np.sin(1.3),np.cos(1.3)]])
    print matrix_has_property(A,'hermitian')    
    print matrix_has_property(A,'square')
    print matrix_has_property(A,'unitary')        
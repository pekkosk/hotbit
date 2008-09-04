import numpy as nu

def separate_symbols(string):
    """ From given string, separate element symbols AuC=>Au,C"""
    if string[1].islower():
        return string[:2],string[2:]
    else:
        return string[0],string[1:]
    
def make_hermitian(H):
    """ Make matrix H hermitian using fast numpy routines. """
    H2=H+H.transpose().conjugate()
    return H2-H.diagonal()*nu.identity(H.shape[0])
        
def make_derivative_antihermitian(dH):
    """ Make vector matrix H antihermitian using fast numpy routines. """
    dH2=nu.zeros_like(dH)
    for a in range(3):
        dH2[:,:,a]=dH[:,:,a]-dH[:,:,a].transpose().conjugate()
    return dH2
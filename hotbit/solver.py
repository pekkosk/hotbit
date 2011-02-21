# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

from scipy.linalg import eig
from scipy.linalg import eigh
import numpy as np
from numpy.linalg import solve
from box.buildmixer import BuildMixer
from weakref import proxy
from random import randint

# Wrapper for the LAPACK dsygvd, zhegvd solvers
from _hotbit import geig

class Solver:
    def __init__(self,calc):
        self.calc=proxy(calc)
        self.maxiter=calc.get('maxiter')
        self.mixer=BuildMixer(calc.mixer)
        self.SCC=calc.get('SCC')
        self.norb=self.calc.el.norb
        self.iterations=None
        self.iter_history=[]

    def __del__(self):
        pass

    def get_nr_iterations(self):
        return self.iterations

    def get_iteration_info(self):
        avg, mx, mn=np.mean(self.iter_history), min(self.iter_history), max(self.iter_history)
        return 'Solved %i times; Iterations: avg %.1f, max %i, min %i' %(len(self.iter_history),avg,mx,mn)


    def get_eigenvalues_and_wavefunctions(self, H0, S, H1=None):
        """
        Solve the generalized eigenvalue problem for a fixed electrostatic
        potential, i.e. a single SCC iteration.
        """
        nk, norb = get_HS_shape(H0, S)

        e  = np.zeros((nk, norb))
        wf = np.zeros((nk, norb, norb), dtype=H0.dtype)

        for ik in range(nk):
            if H1 is not None:
                H = H0[ik] + H1*S[ik]
            else:
                H = H0[ik]
            e[ik], wf[ik] = self.diagonalize(H, S[ik])

        return e, wf


    def get_states(self,calc,dq,H0,S):
        """ Solve the (non)SCC generalized eigenvalue problem. """
        st = calc.st
        es = st.es
        mixer = self.mixer
        mixer.reset()
        H1 = None
        #from box.convergence_plotter import ConvergencePlotter
        #convergence_plotter = ConvergencePlotter(self.calc)
        #convergence_plotter.draw(dq)
        for i in range(self.maxiter):
            # diagonalize for all k-points at once
            if self.SCC:
                H1 = es.construct_h1(dq)
            e, wf = self.get_eigenvalues_and_wavefunctions(H0, S, H1=H1)
            st.update(e,wf)

            # If we don't do SCC, stop here
            if not self.SCC:
                break

            dq_out=st.get_dq()
            done,dq=mixer(dq,dq_out)
            #convergence_plotter.draw(dq)
            if i%10 == 0:
                self.calc.get_output().flush()
            if self.calc.get('verbose_SCC'):
                mixer.echo(self.calc.get_output())
            if done:
                self.iterations=i
                self.iter_history.append(i)
                break
            if i==self.maxiter-1:
                mixer.out_of_iterations(self.calc.get_output())
                #if self.calc.get('verbose_SCC'):
                #    convergence_plotter.show()
                raise RuntimeError('Out of iterations.')
        if self.calc.get('verbose_SCC'):
            mixer.final_echo(self.calc.get_output())
        return st.e,st.wf


    def diagonalize(self,H,S):
        """ Solve the eigenstates. """
        if True:
            # via C wrapper
            self.calc.start_timing('LAPACK eigensolver')
            e, wf = geig(H,S)
            self.calc.stop_timing('LAPACK eigensolver')
            wf = wf.transpose()
            
            # if S happened not to be positive definite, LAPACK
            # results in wrong norm. 
            self.calc.start_timing('Check norm (remove?)')
            norms = ( np.dot(wf.conj(),S)*wf ).sum(axis=1)
            maxdev = np.abs(norms-1).max() 
            if maxdev>1E-6:
                eval,efunc = eigh(S)
                evmin =  eval.min()
                if evmin<0:
                    raise AssertionError('Eigenfunction norm from LAPACK %.4f. Minimum eigenvalue of S is %.4f - overlap matrix is not positive definite.' %(maxdev,evmin))
                else:
                    raise AssertionError('Eigenfunction norm from LAPACK %.4f, but overlap matrix still appears positive definite.' %(maxdev))
            self.calc.stop_timing('Check norm (remove?)')
            
                
        if False:
            #raise NotImplementedError('Not checked for complex stuff')
            # using numpy lapack_lite
            e, wf = eig(H,S)
            e = e.real
            order = e.argsort()
            e = e[order]
            for i in range(self.norb):
                wf[:,i]=wf[order,i]
            for i in range(self.norb): #normalize properly
                wf[i]=wf[i]/np.sqrt( np.abs(np.dot(wf[i],np.dot(S,wf[:].conj()))) )
        
        return e,wf

###

def get_HS_shape(H, S):
    """
    Return the number of k-points and the number of orbitals given a
    Hamiltonian and an overlap matrix.
    """
    nk    = H.shape[0]
    norb  = H.shape[1]
    assert H.shape[2] == norb
    assert S.shape[0] == nk
    assert S.shape[1] == norb
    assert S.shape[2] == norb

    return nk, norb

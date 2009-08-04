# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

from scipy.linalg import eig
import numpy as nu
from hotbit.fortran.eigensolver import geig, geigc
from numpy.linalg import solve
from box.buildmixer import BuildMixer
from weakref import proxy

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
        avg, mx, mn=nu.mean(self.iter_history), min(self.iter_history), max(self.iter_history)
        return 'Solved %i times; Iterations: avg %.1f, max %i, min %i' %(len(self.iter_history),avg,mx,mn)

    #def get_states(self,st,dq,H0,S,count):
        #""" Solve the (non)SCC generalized eigenvalue problem. """
        #self.es=self.calc.es
        #self.norb=len(self.calc.el)
        #if self.SCC:
            #self.es.construct_tables()
        #if True:
            #mixer=AndersonMixer(self.beta,self.mmax,self.limit,chop=0.2)
        #for i in range(self.maxiter):
            #if self.SCC:
                #H=H0+self.es.construct_H1(dq)*S
            #else:
                #H=H0
            #e,wf=self.diagonalize(H,S)
            #st.update(e,wf)
            #if not self.SCC:
                #break
            #dq_out=st.get_dq()
            #done,dq=mixer(dq,dq_out)

            #if self.verbose:
                #mixer.echo()
            #if done:
                #self.iterations=i
                #self.iter_history.append(i)
                #break
            #if i==self.maxiter-1:
                #mixer.out_of_iterations(self.calc.get_output())
                #raise RuntimeError('Out of iterations.')
        #if self.verbose:
            #mixer.final_echo()
        #return st.e,st.wf


    def get_states(self,calc,dq,H0,S,count):
        """ Solve the (non)SCC generalized eigenvalue problem. """
        st = calc.st
        es = st.es
        mixer = self.mixer
        mixer.reset()
        #from box.convergence_plotter import ConvergencePlotter
        #convergence_plotter = ConvergencePlotter(self.calc)
        #convergence_plotter.draw(dq)
        e=nu.zeros((st.nk,self.norb))
        wf=nu.zeros((st.nk,self.norb,self.norb),complex)
        for i in range(self.maxiter):
            # diagonalize for all k-points at once
            for ik in range(st.nk):
                if self.SCC:
                    H=H0[ik,:,:] + es.construct_h1(dq)*S[ik,:,:]
                else:
                    H=H0[ik,:,:]
                e[ik,:], wf[ik,:,:] = self.diagonalize(H,S[ik,:,:])
                self.H, self.S=H,S
            st.update(e,wf)
            if not self.SCC:
                break
            
            #raise NotImplementedError
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
        self.calc.start_timing('eigensolver')
        if True:
            # via Fortran wrapper
#            e,wf=geig(self.H,self.S,self.norb)
            e, wf = geigc(H,S,self.norb)
        if False:
            raise NotImplementedError('Not checked for complex stuff')
            # using numpy lapack_lite
            e,wf=eig(self.H,self.S)
            e=e.real
            order=e.argsort()
            e=e[order]
            for i in range(self.norb):
                wf[i,:]=wf[i,order]
            for i in range(self.norb): #normalize properly
                wf[:,i]=wf[:,i]/nu.sqrt( nu.dot(wf[:,i],nu.dot(self.S0,wf[:,i])) )
        self.calc.stop_timing('eigensolver')
        return e,wf







import numpy as nu
from box.dummymixer import DummyMixer

class PulayMixer(DummyMixer):
    """ A class for Pulay mixing, as described in
        https://wiki.fysik.dtu.dk/gpaw/documentation/densitymix/densitymix.html#densitymix """


    def __init__(self, mixing_constant=0.3, memory=8, convergence=1e-3):
        DummyMixer.__init__(self, mixing_constant, convergence)
        self.name = 'Pulay'
        self.memory = memory
        self.initialized = False
        self.A = nu.zeros((self.memory, self.memory))
        self.alfa = nu.zeros(self.memory)


    def _initialize(self, xi):
        """ Allocate fixed-sized arrays. """
        self.R = nu.zeros((self.memory, len(xi)))
        self.rho = nu.zeros((self.memory, len(xi)))
        self.initialized = True


    def  __call__(self, xi, yi):
        simple = False
        self.it += 1
        if not self.initialized:
            self._initialize(xi)
        lim = min(self.it, self.memory)
        self.rho[1:lim] = self.rho[0:lim-1]
        self.rho[0] = xi
        self.R[1:lim] = self.R[0:lim-1]
        self.R[0] = yi - xi
        for i in range(lim):
            for j in range(lim):
                self.A[i,j] = nu.dot(self.R[j], self.R[i])
        try:
            inv_A = nu.linalg.inv(self.A[:lim,:lim])
        except:
            # problem in taking the inverse of the matrix A,
            # use simple mixing in this step.
            simple = True
        if simple:
            xb = (1-self.beta)*xi + self.beta*yi
        else:
            for i in range(lim):
                up = nu.sum(inv_A[:lim, i])
                down = nu.sum(inv_A[:lim,:lim])
                self.alfa[i] = up/down
            xb = nu.zeros_like(xi)
            for i in range(lim):
                xb += self.alfa[i]*(self.rho[i] + self.beta*self.R[i])
        fmax = max(abs(self.R[0]))
        self.fmax.append(fmax)
        if fmax < self.convergence:
            return True, xb
        else:
            return False, xb


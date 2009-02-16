import numpy as nu

class PulayMixer:
    """ A class for Pulay mixing, as described in
        https://wiki.fysik.dtu.dk/gpaw/documentation/densitymix/densitymix.html#densitymix """

    def __init__(self, beta, mmax, limit):
        raise NotImplementedError("This mixer has not been tested well. Comment this line in order to use it.")
        self.beta = beta
        self.mmax = mmax
        self.limit = limit
        self.initialized = False
        self.it = 0

        self.A = nu.zeros((mmax, mmax))
        self.alfa = nu.zeros(mmax)

    def  __call__(self, xi, yi):
        simple = True
        self.it += 1
        if not self.initialized:
            self._initialize(xi)
        M = self.mmax
        self.rho[1:M] = self.rho[0:M-1]
        self.R[1:M] = self.R[0:M-1]
        self.alfa[1:M] = self.alfa[0:M-1]
        self.rho[0] = xi
        self.R[0] = yi - xi
        lim = min(self.it, M)
        for i in range(lim):
            for j in range(lim):
                # FIXME no need to calculate everything again
                self.A[i,j] = nu.dot(self.R[j], self.R[i])
        try:
            inv_A = nu.linalg.inv(self.A[:lim,:lim])
            simple = False
        except:
            return False, (1-self.beta)*xi + self.beta*yi
        for i in range(lim):
            up = nu.sum(inv_A[:lim, i])
            down = nu.sum(inv_A[:lim,:lim])
            self.alfa[i] = up/down
        xb = nu.zeros_like(xi)
        R_opt = nu.zeros_like(self.R[0])
        for i in range(lim):
            xb += self.alfa[i]*(self.rho[i] + self.beta*self.R[i])
            R_opt += self.alfa[i]*self.R[i]
        Rmax = max(abs(self.R[0]))
        self.Rmax.append(Rmax)
        #if self.it < 2:
        #    return False, (1-self.beta)*xi + self.beta*yi
        #elif 2 <= self.it <= M:
        if self.it <= M:
            return False, xb
        else:
            if Rmax < self.limit:
                return True, xb
            else:
                return False, xb

    def _initialize(self, xi):
        self.R = nu.zeros((self.mmax, len(xi)))
        self.Rmax = []
        self.rho = nu.zeros((self.mmax, len(xi)))
        self.initialized = True

    def echo(self):
        return "Pulay: iter: %i   fmax:  %0.12f" % (self.it, self.Rmax[-1])

    def final_echo(self):
        return self.it

    def out_of_iterations(self, out):
        print >> out, "Pulay mixer out of iterations!"
        for it, fmax in enumerate(self.Rmax):
            print >> out, it, fmax

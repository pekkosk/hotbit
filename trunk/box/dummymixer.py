class DummyMixer:
    """ A dummy mixer that all mixer class should inherit. """


    def __init__(self, beta, convergence):
        self.name = 'DummyMixer'
        self.beta = beta
        self.convergence = convergence
        self.fmax = []
        self.it = 0


    def __call__(self, xi, yi):
        return False, yi


    def get(self, key):
        return self.__dict__[key]


    def echo(self, out):
        """ Say something about the progress of iteration. """
        print >> out, "%s: iter %i  fmax %0.12f" %  (self.name, self.it, self.fmax[-1])


    def final_echo(self, out):
        print >> out, self.it


    def out_of_iterations(self, out):
        """ Print out-of-iterations info. """
        print>>out, '%s mixer out of iterations!' % self.name
        print>>out, 'History of maximum residuals:'
        for it,fmax in enumerate(self.fmax):
            print>>out, it,fmax


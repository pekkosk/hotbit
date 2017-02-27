from __future__ import print_function

class DummyMixer:
    """ A dummy mixer that all mixer class should inherit. """


    def __init__(self, beta, convergence):
        self.name = 'DummyMixer'
        self.beta = beta
        self.convergence = convergence
        self.reset()


    def reset(self):
        self.fmax = []
        self.it = 0


    def __call__(self, xi, yi):
        return False, yi


    def get(self, key):
        return self.__dict__[key]


    def echo(self, out):
        """ Say something about the progress of iteration. """
        print("%s: iter %i  fmax %0.12f" %  (self.name, self.it, self.fmax[-1]), file=out)
        out.flush()


    def final_echo(self, out):
        print(self.it, file=out)


    def out_of_iterations(self, out):
        """ Print out-of-iterations info. """
        print('%s mixer out of iterations!' % self.name, file=out)
        print('History of maximum residuals:', file=out)
        for it,fmax in enumerate(self.fmax):
            print(it,fmax, file=out)


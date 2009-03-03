import numpy as nu

class ConvergencePlotter:
    """ Show Mulliken excess charges generated during the iterative
    solving cycle. """

    def __init__(self, calc):
        M = 10 # how many steps to show
        N = len(calc.el)
        self.maxiter = min(calc.st.solver.maxiter, M)
        self.title = calc.el.get_name() + ", %i last excess Mulliken populations, (red=latest, green=oldest)" % self.maxiter
        self.N1 = self.maxiter/2
        if self.maxiter % 2 == 0:
            self.N2 = self.maxiter/2
        else:
            self.N2 = self.maxiter/2+1
        # colors from green to yellow
        colors1 = [(n,1,0) for n in nu.linspace(0,1,self.N1)]
        # colors from yellow to red
        colors2 = [(1,1-n,0) for n in nu.linspace(0,1,self.N2)]
        self.colors = colors1 + colors2
        self.it = 0
        self.points = nu.zeros((self.maxiter,N))

    def color(self):
        """ Return next color. """
        ret = self.colors[self.it]
        self.it += 1
        return tuple(ret)

    def draw(self, dq):
        """ Add points to be shown later. """
        self.points[1:self.maxiter,:] = self.points[0:self.maxiter-1,:]
        self.points[0] = dq

    def show(self):
        """ Show the last M points that are still in the memory. """
        try:
            import pylab
        except ImportError:
            print "Could not import pylab, cannot print the Mulliken excess charges."
            return
        for p in self.points:
            c = self.color()
            pylab.scatter(range(len(p)), p, color=c)
        pylab.title(self.title)
        pylab.xlabel('atom index')
        pylab.ylabel('Mulliken excess population')
        pylab.savefig('mulliken_charge_fluctuations_%i.png' % self.maxiter)
        pylab.show()

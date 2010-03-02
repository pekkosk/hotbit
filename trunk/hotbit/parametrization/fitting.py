from hotbit import Element
from hotbit import Calculator
import numpy as nu
from ase.units import Bohr,Hartree
from box import mix
from box import Atoms
from box.interpolation import Function
from ase import read, PickleTrajectory
import pylab as pl
from copy import copy
from sys import stdout





class RepulsiveFitting:

    def __init__(self,symbol1,symbol2,r_cut,s=None,k=3,txt=None,tol=1E-5): 
        """
        Class for fitting the short-range repulsive potential.
        
        
        Fitting uses eV and Angstrom also internally, only the
        output file (.par) is in Hartrees and Bohrs  
        
        Parameters:
        ===========
        symbol1:        chemical symbol for the first element
        symbol2:        chemical symbol for the second element
        txt:            output filename or None for stdout
        r_cut:          the repulsion cutoff 
        s:              smoothing parameter. If None, use s = N - nu.sqrt(2*N)
                        where N is the number of data points.
        k:              order of spline, cubic by default.
                        Uses smaller order if not enough points to fit V_rep'(R)                      
        tol:            tolerance for distances still considered the same
        
        
        Usage:
        ======        
        1. Initialize class
           * rep = RepulsiveFitting('Au','Au',r_cut=3.3,s=100)
    
        2. Collect data from structures. The data collected is points r_i and V_rep'(r_i),
           that is, repulsive distance and the force. Use the methods 'append_*'.
           * rep.append_dimer(weight=0.5,calc=calc0,R=2.49,comment='Au2')
           * rep.append_energy_curve(weight=1.0,calc=calc0,traj='dimer_curve.traj',label='DFT dimer',comment='dimer curve')
    
        3. Given the set of points [r_i,V_rep'(r_i)], fit a spline with given order.
           * fit()
           Fitting will produce a spline-interpolated V_rep'(r), which is then integrated
           to given spline-interpolated V_rep(r).
    
        4. Output repulsion into a file and plot the repulsion
           * rep.write_par('Au_Au_no_repulsion.par',filename='Au_Au_repulsion.par')
           * rep.plot('AuAu_repulsion.pdf')
    
        """
        self.elm1=Element(symbol1)
        self.elm2=Element(symbol2)
        self.sym1=symbol1
        self.sym2=symbol2
        self.r_cut = r_cut
        self.s = s
        self.k = k
        self.tol = tol
        
        self.r_dimer = None
        self.deriv=[]
        self.comments=''
        self.scale=1.025                    # scaling factor for scalable systems
        self.structures = []
        self.v=None

        if txt==None:
            self.txt=stdout
        else:
            self.txt=open(txt,'a')
        print>>self.txt, 'Fitting repulsion curve between %s and %s' % (self.sym1, self.sym2)           



    def __call__(self,r,der=0):
        """ 
        Return repulsion or its derivative.
        
        parameters:
        ===========
        r:            radius
        der:          der=0 for V_rep(r), der=1 for V_rep'(r)
        """
        return self.v(r,der=der)


    def plot(self, filename=None):
        """ 
        Plot vrep and derivative together with fit info. 
        
        parameters:
        ===========
        filename:     graphics output file name
        """
        r=nu.linspace(0,self.r_cut)
        v=[self(x,der=0) for x in r]
        vp=[self(x,der=1) for x in r]
        rmin=0.9*min([d[0] for d in self.deriv])

        fig=pl.figure()
        pl.subplots_adjust(wspace=0.25)
        
        # Vrep
        pl.subplot(1,2,1)
        pl.ylabel(r'$V_{rep}(r)$  (eV)')
        pl.xlabel(r'$r$  ($\AA$)')
        pl.axvline(x=self.r_cut,c='r',ls=':')
        pl.plot(r,v)
        pl.ylim(ymin=0,ymax=self(self.r_cut/2))
        pl.xlim(xmin=self.r_cut/2, xmax=self.r_cut)

        # Vrep'
        pl.subplot(1,2,2)
        pl.ylabel(r'$dV_{rep}(r)/dr$ (eV/$\AA$)')
        pl.xlabel(r'$r$ ($\AA$)')
        pl.plot(r,vp,label=r'$dV_{rep}(r)/dr$')
        for s in self.deriv:
            pl.scatter( [s[0]],[s[1]],s=100*s[2],c=s[3],label=s[4])
        pl.axvline(x=self.r_cut,c='r',ls=':')
        if self.r_dimer!=None:
            pl.axvline(x=self.r_dimer,c='r',ls=':')
        xmin = 0.8*self.r_dimer
        xmax = 1.2*self.r_cut
        ymin = min([ point[1] for point in self.deriv ])
        ymax = nu.abs(ymin)*0.2
        pl.axhline(0,ls='--',c='k')
        pl.text(self.r_dimer, ymax, r'$r_{dimer}$')
        pl.text(self.r_cut, ymax, r'$r_{cut}$')
        pl.xlim(xmin=xmin, xmax=xmax)
        pl.ylim(ymin=ymin, ymax=ymax)
        #pl.subtitle('Fitting for %s and %s' % (self.sym1, self.sym2))
        pl.rc('font', size=10)
        pl.legend(loc=4)
        file = '%s_%s_repulsion.pdf' % (self.sym1, self.sym2)
        if filename!=None:
            file=filename
        pl.savefig(file)
        pl.clf()


    def add_comment(self,s=None):
        """ 
        Append some comment for par-file.
        
        These comments will end up in Hotbit's output logfile each time
        fitted repulsion is used in calculations. For this reason,
        use as short and concise comments as possible. 
        
        parameters:
        ===========
        s:            comment as a string
        """
        if s in [None, '']:
            return
        add='|'
        if len(self.comments)==0:
            add=''
        self.comments+=add+s


    def fit(self): 
        """
        Fit spline into {r, V_rep'(r)} -data.
        
        Use cubic spline by default
        
        parameters:
        ===========
 
        """
        from scipy.interpolate import splrep, splev
        self.k = min(len(self.deriv)-1+1,self.k)
        
        x = nu.array([self.deriv[i][0] for i in range(len(self.deriv))])
        y = nu.array([self.deriv[i][1] for i in range(len(self.deriv))])
        w = nu.array([self.deriv[i][2] for i in range(len(self.deriv))])
        # sort values so that x is in ascending order
        indices = x.argsort()
        x = x[indices]
        y = y[indices]
        w = w[indices]
        x, y, w = self._group_closeby_points(x,y,w)
        # use only points that are closer than r_cut
        indices = nu.where(x < self.r_cut)
        x = list(x[indices])
        y = list(y[indices])
        w = list(w[indices])
        # force the spline curve to go to zero at x=r_cut
        x.append(self.r_cut)
        y.append(0)
        w.append(1e3*max(w))
        if self.s == None:
            # from documentation of splrep in scipy.interpolate.fitpack
            self.s = len(x) - nu.sqrt(2*len(x))
            
        print>>self.txt, "\n  Fitting smoothing spline with parameters"
        print>>self.txt, "  k=%i, s=%0.4f, r_cut=%0.4f\n" %(self.k, self.s, self.r_cut)
        tck = splrep(x, y, w, s=self.s, k=self.k)

        def dv_rep(r):
            return splev(r, tck)

        v_rep = self._integrate_vrep(dv_rep, self.r_cut)

        def potential(r, der=0):
            if der == 0:
                return v_rep(r)
            elif der == 1:
                return dv_rep(r)
            else:
                raise NotImplementedError("Only 0th and 1st derivatives")
        self.v = potential
        

    def _group_closeby_points(self, x, y, w):
        """
        If there are many y-values with almost the same x-values,
        it is impossible to make spline fit to these points.
        For these points the y will be the weighted average of
        the y-points and the weight is the sum of the weights of
        averaged y-points.
        """
        accuracy = 4 # the number of decimals to maintain
        pseudo_x = nu.array(x*10**accuracy, dtype=int)
        groups = nu.zeros(len(x), dtype=int)
        g = 0
        for i in range(1,len(pseudo_x)):
            if pseudo_x[i] != pseudo_x[i-1]:
                groups[i] = groups[i-1] + 1
            else:
                groups[i] = groups[i-1]
        new_x = []
        new_y = []
        new_w = []
        for g in range(max(groups)+1):
            same = nu.where(groups == g)
            new_x.append(nu.average(x[same]))
            new_y.append(nu.dot(y[same],w[same])/nu.sum(w[same]))
            new_w.append(nu.sum(w[same]))
        return nu.array(new_x), nu.array(new_y), nu.array(new_w)


    def write_par(self, inputpar, filename=None):
        """
        Write the full par-file to file.
        
        parameters:
        ===========
        inputpar:     the par-file where the repulsion is appended
        filename:     output file
        """
        from time import asctime
        import shutil

        if filename==None:
            filename = 'repulsion_'+inputpar
            
        shutil.copy(inputpar, filename)
        f = open(filename, 'a')
        # add comments
        print >> f, "repulsion_comment="
        print >> f, "%s\nparameters r_cut = %0.4f Ang, s = %0.4f, k = %3i" % (asctime(),self.r_cut, self.s, self.k)
        if len(self.structures)>1:
            print >> f, "The systems used to produce this fit:"
            for data in self.structures:
                print >> f, "%20s %3s" % (data['filename'], data['charge'])
        if len(self.comments) > 0:
            print >> f, self.comments            
        print >> f, '\n\nrepulsion='
        for r in nu.linspace(0.1, self.r_cut, 100):
            print >> f, r/Bohr, self(r)/Hartree
        f.close()


    def _integrate_vrep(self, dv_rep, r_cut, N=100):
        """
        Integrate V'_rep(r) from r_cut to zero to get the V_rep(r)
        """
        from box.interpolation import SplineFunction
        from scipy.integrate import quadrature
        r_g = nu.linspace(r_cut, 0, N)
        dr = r_g[1] - r_g[0]
        v_rep = nu.zeros(N)
        for i in range(1,len(r_g)):
            v_rep[i] = v_rep[i-1]
            val, err = quadrature(dv_rep, r_g[i-1], r_g[i], tol=1.0e-12, maxiter=50)
            v_rep[i] += val
        # SplineFunction wants the x-values in ascending order
        return SplineFunction(r_g[::-1], v_rep[::-1])


    def _set_calc(self,atoms,calc):
        """
        Set calculator for given atoms.
        """
        if type(atoms)==type(''):
            a = read(atoms)
        else:
            a = atoms.copy()
        c = copy(calc)
        c.set('txt','-')
        a.set_calculator(c)
        return a,c

    #
    #       Fitting methods
    #             
    def append_point(self,weight,R,dvrep,comment=None,label='point',color='g'):
        """ Add point to vrep'-fitting.
        
        parameters:
        ===========
        weight:    relative weight (inverse of deviation)
        R:         radius, in Angstroms
        vder:      value of the derivative at R
        color:     pylab color for plotting
        info:      label in plotting
        """
        self.deriv.append([R,dvrep,weight,color,label])
        self.add_comment(comment)


    def append_scalable_system(self,weight,calc,atoms,comment=None,label=None,color='m'):
        """ 
        Use scalable equilibrium system in repulsion fitting. 
        
        parameters:
        ===========
        weight:        fitting weight 
        calc:          Hotbit calculator (remember charge and k-points)
        atoms:         filename or ase.Atoms
        comment:       fitting comment for par file
        label:         plotting label
        color:         plotting color
        """
        atoms, calc = self._set_calc(atoms,calc)
        if label is None:
            label=atoms.get_name()+'_scalable'
        
        e1 = atoms.get_potential_energy()
        R, N = self._get_repulsion_distances(calc)
        atoms.set_cell( atoms.get_cell()*self.scale, scale_atoms=True )
        e2 = atoms.get_potential_energy()

        dEdr=(e2-e1)/(self.scale*R-R)
        self.append_point(weight,R,-dEdr/N,comment,label,color)
        print>>self.txt, '\nAdding a scalable system %s with %i bonds at R=%.4f.' %(atoms.get_name(),N,R)


    def _get_repulsion_distances(self,calc):
        """
        Return distances below r_cut for given system in calculator.
        
        return:
        =======
        R:     the mean repulsion distance
        N:     number of bonds
        """
        distances = calc.rep.get_repulsion_distances(self.sym1,self.sym2,self.r_cut)
        R = distances.mean()
        if len(distances)==0:
            return 0.0,distances
        if distances.max()-distances.min() > self.tol:
            atoms = calc.get_atoms()
            raise AssertionError('Bond lengths in %s are not the same' %atoms.get_name() )
        N = len(distances)
        return R,N        


    def append_dimer(self,weight,calc,R,comment=None,label='dimer',color='r'):
        """ 
        Use dimer bond length in fitting. 
        
        parameters:
        ===========
        weight:    relative weight (inverse of deviation)
        calc:      Hotbit calculator used in calculation (remember Gamma-point and charge)
        R:         dimer bond length (Angstroms)
        comment:   fitting comment for par-file
        label:     plotting label
        color:     plotting color        
        """
        self.r_dimer = R
        atoms = Atoms([self.sym1,self.sym2],[(0,0,0),(R,0,0)],pbc=False)
        atoms.center(vacuum=5)
        self.append_scalable_system(weight,calc,atoms,comment=comment,label=label,color=color)


    def append_energy_curve(self,weight,calc,traj,comment=None,label='energy curve', color='y'):
        """
        Calculates the V'rep(r) from a given ase-trajectory.
        
        The trajectory can be anything, as long as the ONLY missing energy
        from DFTB calculation is N*V_rep(R). Hence

                     E_DFT'(R) - E_wr'(R)
        V_rep'(R) =  ------------------ ,
                             N

        where R is the nn. distance,N is the number of A-B pairs taken into account,
        and E_wr(R) = E_bs(R) + E_coul(R) is the DFTB energy without repulsion.
        At least 3 points in energy curve needed, preferably more.

        parameters:
        ===========        
        weight:              Fitting weight (inverse of deviation)
        calc:                Hotbit calculator (remember charge and k-points)
        traj:                ASE trajectory for energy curve
        comment:             fitting comment for par-file        
        label:               Name that is used in the fitting plot.
        color:               Color that is used in the fitting plot.
        """
        from box.interpolation import SplineFunction

        print>>self.txt, "\nAppending energy curve data from %s..." %traj
        traj = PickleTrajectory(traj)
        Edft, Ewr, N, R = [], [], [], []
        if len(traj)<3:
            raise AssertionError('At least 3 points in energy curve required.')
        for a in traj:
            atoms, c = self._set_calc(a,calc)
            e = atoms.get_potential_energy()
            r, n = self._get_repulsion_distances(c)
            if n>0:
                Edft.append( a.get_potential_energy() )
                Ewr.append( e )
                R.append(r)
                N.append(n)
        Edft = nu.array(Edft)
        Ewr = nu.array(Ewr)
        N = nu.array(N)
        R = nu.array(R)
            
        # sort radii because of spline
        indices = R.argsort()
        R    = R[indices]
        Edft = Edft[indices]
        Ewr  = Ewr[indices]
        vrep = SplineFunction(R, (Edft-Ewr)/N, k=3, s=0)

        for i, r in enumerate(R):
            if i > 0: 
                label='_nolegend_'
                com = None
            else:
                com = comment
            self.append_point(weight/nu.sqrt(len(R)),r, vrep(r,der=1), com, label, color)
        print>>self.txt, "Appended %i points around R=%.4f...%.4f" %(len(N),R.min(),R.max())


    def append_homogeneous_cluster(self,weight,calc,atoms,comment=None,label='cluster',color='b'):
        """
        Use homonuclear cluster in fitting, even having different bond lengths.
        
        Construct repulsive forces so that residual forces |F_DFT - (F_WR+F_rep)|,
        where F_DFT are DFT forces (zero if cluster in equilibrium), F_WR are
        DFTB forces without repulsion, and F_rep are the repulsive forces. That is,
        minimize the function
        
           sum_i |F_DFT_i - F_WR_i - F_rep_i|^2
           
        with respect a and b, where V_rep'(R) = a + b*(r-r_cut). Then, add fitting points
        from rmin to rmax, where these values span all pair distances below r_cut
        within the cluster.
         
        Only finite, non-periodic systems can be used.
        
        parameters:
        ===========
        weight:        relative weight
        calc:          Hotbit calculator (remember charge and no k-points)
        atoms:         filename or ASE.Atoms instance
        comment:       fitting comment for par-file
        label:         plotting label
        color:         plotting color
        """
        if type(atoms)==type(''):
            atoms = read(atoms)

        N = len(atoms)
        try:
            f_DFT = atoms.get_forces()
            print>>self.txt, "    Use forces"
        except:
            f_DFT = nu.zeros((N,3))
            print>>self.txt, "    No forces (equilibrium cluster)"
            
        atoms, calc = self._set_calc(atoms,calc)
        print>>self.txt, "\nAppending homogeneous cluster %s..." % atoms.get_name()
        
        f_WR = atoms.get_forces()
        distances = calc.rep.get_repulsion_distances(self.sym1,self.sym2,self.r_cut)
        rmin, rmax = distances.min(), distances.max()
        
        def dvrep(r,p):
            """ Auxiliary first-order polynomial for repulsion derivative """
            return p[0]+p[1]*(r-self.r_cut)
            
        
        def to_minimize(p,atoms,fdft,fwr):
            """ Function sum_I |F_DFT_I - F_TB_I|^2 to minimize. """
            N = len(atoms)
            pos = atoms.get_positions()
            resid = nu.zeros((N,3))
            frep  = nu.zeros((N,3))
            for i in range(N):
                for j in range(N):
                    if i==j: continue
                    vec = pos[j]-pos[i]
                    d = nu.linalg.norm(vec)
                    if d>self.r_cut: 
                        continue
                    else:
                        frep[i] += dvrep(d,p)*vec/d
            resid = fdft - ( fwr + frep )
            return sum([ nu.linalg.norm(resid[i]) for i in range(N) ])                    
        
        from scipy.optimize import fmin
        p = fmin( to_minimize,[-1.0,5.0],args=(atoms,f_DFT,f_WR),xtol=1E-5,ftol=1E-5 )
        print>>self.txt, '   Cluster: V_rep(R)=%.6f + %.6f (r-%.2f)' %(p[0],p[1],self.r_cut)
      
        np = 6
        rlist = nu.linspace(rmin,rmax,np)
        for i,r in enumerate(rlist):
            if i==0:
                com = comment
            else:
                label = '_nolegend_'
                com = None
            self.append_point(weight/nu.sqrt(np), r, dvrep(r,p), com, label, color)


    def write_fitting_data(self, filename):
        import pickle
        f = open(filename,'w')
        pickle.dump(self.deriv, f)
        pickle.dump(self.structures, f)
        pickle.dump(self.comments, f)
        f.close()


    def load_fitting_data(self, filename):
        import pickle
        f = open(filename,'r')
        self.deriv = pickle.load(f)
        self.structures = pickle.load(f)
        self.comments = pickle.load(f)
        f.close()


    def _get_trajs_for_fitting(self):
        return self.structures
    
    

class ParametrizationTest:
    """
    A tool to examine how well your parametrization agrees with
    given ase-trajectories.

    trajectories:  list of trajectories you want to compare
    charges:       the charges of the systems in trajectories
    """

    def __init__(self, rf, pars):
        from copy import copy
        from hotbit import Hotbit
        self.pars = pars
        self.trajectories = []
        self.calculators = []
        for data in rf._get_trajs_for_fitting():
            filename = data['filename']
            del data['filename']
            c = Hotbit()
            for key, value in data.iteritems():
                c.__dict__[key] = data[key]
            self.trajectories.append(filename)
            self.calculators.append(c)
        self.points = []
        self.ref_points = []
        self.colors = ['cyan','red','orange','#8DEE1E','magenta','green','black']

    def norm_to_isolated_atoms(self, atoms):
        """
        Return the constant that can be used to calculate
        the binding energy of the system.
        """
        delta_E = 0
        for atom in atoms:
            delta_E -= self.E_free[atom.symbol]
        return delta_E


    def get_isolated_energies(self, trajs, par):
        """
        Return the energies of an isolated atoms.
        """
        elements = []
        energies = {}
        for t in trajs:
            traj = PickleTrajectory(t)
            for atom in traj[0]:
                if not atom.symbol in elements:
                    elements.append(atom.symbol)
        el1, el2 = par.split("_")[0:2]
        for el in elements:
            ss = "%s%s" % (el, el)
            if el1 == el2 and el1 == el:
                tables = {ss:par, 'others':'default'}
                calc = Calculator(SCC=True, tables=tables)
            else:
                calc = Calculator(SCC=True)
            atoms = Atoms(ss, ((0,0,0),(200,0,0)))
            atoms.center(vacuum=100)
            atoms.set_calculator(calc)
            energies[el] = atoms.get_potential_energy() / 2
        return energies


    def compare(self):
        """
        Make a comparison for all the systems.
        """
        for i_par in range(len(self.pars)):
            self.compare_with_par(i_par)


    def compare_with_par(self, i_par):
        """
        Make a comparison to all trajectories with given parameter-file.
        The i_par is the index to the self.pars.
        """
        par = self.pars[i_par]
        self.E_free = self.get_isolated_energies(self.trajectories, par)
        temp = par.split('_')
        symbols = "%s%s" % (temp[0],temp[1])
        tables = {symbols:par, 'others':'default'}
        for i_traj, calc in zip(range(len(self.trajectories)), self.calculators):
            pl.figure(i_traj)
            pl.title(self.trajectories[i_traj])
            if i_par == 0:
                self.plot_ref(i_traj)
            self.compare_trajectory(i_traj, calc, tables, i_par)


    def compare_trajectory(self, i_traj, calc, tables, i_par):
        """
        Calculate the energies for the frames in the trajectory
        and plot them.
        """
        frames = []
        energies = []
        trajectory = PickleTrajectory(self.trajectories[i_traj])
        for i, image in enumerate(trajectory):
            e_tb = None
            try:
                atoms = Atoms(image)
                c = copy(calc)
                c.tables = tables
                atoms.set_calculator(c)
                e_tb = atoms.get_potential_energy()
            except Exception, ex:
                print>>self.txt, ex
            if e_tb != None:
                energies.append(e_tb)
                frames.append(i)
        delta_E = self.norm_to_isolated_atoms(trajectory[0])
        for i in range(len(energies)):
            energies[i] += delta_E
        self.plot(frames, energies, i_traj, tables, i_par)


    def plot_ref(self, i_traj):
        """
        Plot the energies of a given trajectory as a function
        of the frame number.
        """
        e_dft = []
        traj = PickleTrajectory(self.trajectories[i_traj])
        for image in traj:
            e_dft.append(image.get_total_energy())
        pl.plot(e_dft, c='blue', label='DFT-energies')


    def plot(self, frames, points, i_traj, tables, i_par):
        par = self.pars[i_par]
        color = self.colors[i_par % len(self.colors)]
        pl.plot(frames, points, c=color, label='TB-%s' % par)
        pl.xlabel('frame #')
        pl.ylabel('Energy (eV)')
        pl.legend()


    def run(self):
        """
        Make all the comparisons with given trajectories and parameter
        files and show the results.
        """
        self.compare()
        pl.show()


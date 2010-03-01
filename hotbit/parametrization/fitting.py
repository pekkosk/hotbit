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
        R:         dimer bond length in Angstroms. 
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
        vrep = SplineFunction(R, (Edft-Ewr)/N)

        for i, r in enumerate(R):
            if i > 0: 
                label='_nolegend_'
                com = None
            else:
                com = comment
            self.append_point(weight/nu.sqrt(len(R)),r, vrep(r,der=1), com, label, color)
        print>>self.txt, "Appended %i points around R=%.4f...%.4f" %(len(N),R.min(),R.max())


    def append_homogeneous_cluster(self,weight,calc,atoms,tol=0.005,comment=None,label='homogeneous structure',color=None):
        """
        Use homonuclear cluster in fitting, even having different bond lengths.
        
        Construct repulsive forces so that residual forces are minimized 
        (F_i = \sum_j(dV(|r_ij|)/dR)). System can be stable (no forces
        or forces==0), but it can also be any structure with given forces.
        Only finite, non-periodic systems.
        
        parameters:
        ===========
        weight:        relative weight
        calc:          Hotbit calculator (remember charge and k-points)
        atoms:         filename or ASE.Atoms instance
        tol:           acceptable tolerance in bond lengths (still considered equal)
        comment:       fitting comment for par-file
        label:         plotting label
        color:         plotting color
        """
        atoms = self._set_calc(atoms,calc)
        N = len(atoms)
        if any(atoms.get_pbc()):
            raise AssertionError('Cluster should not be periodic')
        
        print>>self.txt, "\nAppending homogeneous cluster %s..." % atoms.get_name()
               
        epsilon, distances, mask = self._get_matrices(atoms, tol)
        if len(distances) == 0:
            raise RuntimeError("There are no bonds under given cut radius in frame %i." & fr)
        r_hat = self._construct_rhat_matrix(atoms)

        # if the given structure contains forces, use them, otherwise
        # assume that the structure is optimized
        try:
            dft_forces = atoms.get_forces()
            print>>self.txt, "    Use forces"
        except:
            dft_forces = nu.zeros((N,3))
            print>>self.txt, "    No forces (equilibrium cluster)" 
        
        atoms.get_potential_energy()
        residual = dft_forces - atoms.get_forces()

        # use one less point in an array that is given for
        # the minimizer to reduce the number of degrees of freedom
        v_rep_points = nu.zeros(len(distances)-1)

        def residual_forces(v_rep_points):
            # the function that is minimized, returns the sum of
            # the norm of forces acting on each atoms
            res = 0.
            for i in range(N):
                indices = nu.array(mask[i,:]*epsilon[i,:])
                indices -= 1
                ppoints  = v_rep_points[indices]
                f = residual[i].copy()
                f -= nu.dot((mask[i,:]*ppoints),r_hat[i,:])
                res += nu.dot(f,f)
            return res

        print>>self.txt, "    Found %i different bond lengths." % len(v_rep_points)
        v_rep_points, last_res_forces, minimized = self._find_forces(residual_forces, v_rep_points)
        print>>self.txt, "    The sum of the squared norms of the net forces: %0.4f (eV/ang)**2" % (last_res_forces)
        # finally add the missing component
        v_rep_points = list(v_rep_points)
        v_rep_points.insert(0,0)
        if minimized:
            points = []
            for i in range(0,N-1):
                for j in range(i+1,N):
                    if mask[i,j] > 0:
                        # could add a degeneracy factor since there
                        # may be different number of different bonds
                        points.append([distances[epsilon[i,j]], v_rep_points[epsilon[i,j]], last_res_forces])
        else:
            print>>self.txt, "    The minimization of forces did not converge!"
                
        if len(points) > 0:
            points = nu.array(points)
            sigmas = points[:,2] + nu.min(points[:,2])*0.001 # add small value to prevent 1/sigma to go infinity
            inv_sigmas = 1./sigmas
            norm_factor = nu.sqrt(weight**2 / nu.dot(inv_sigmas, inv_sigmas))
            points[:,2] /= norm_factor
            for data in points:
                self.append_point(1/data[2], data[0], data[1], comment, label, color)
                label = '_nolegend_'
            self.add_comment(comment)
                    
                        
                
    def _construct_rhat_matrix(self, structure):
        N = len(structure)
        r_hat = nu.zeros((N,N,3))
        for i in range(N):
            for j in range(N):
                if not i == j:
                    vec = structure.positions[j] -structure.positions[i]
                    norm = nu.linalg.norm(vec)
                    r_hat[i,j] = vec/norm
        return r_hat


    def _get_matrices(self, atoms, tol):
        """
        Construct epsilon matrix that maps the indices (i,j) to a
        single list of distances. If there are many bonds with
        almost same lengths, treat these bonds as there was only
        one of them in order to reduce the degree of freedom
        in the minimization. If the bond is longer that given
        cut radius, that bond is ignored.
        
        parameters:
        ===========
        atoms:        ASE.Atoms instance
        tol:          tolerance for bond lengths
        """
        N = len(atoms)
        distances = []
        index_list = []
        for i in range(N):
            for j in range(N):
                vec = atoms.positions[j] - atoms.positions[i]
                distances.append(nu.linalg.norm(vec))
                index_list.append([i,j])
        distances = nu.array(distances)
        index_list = nu.array(index_list)

        indices = distances.argsort()
        distances = distances[indices]
        index_list = index_list[indices]
        groups = nu.zeros(len(distances), dtype=int)
        group = 0
        for i, d in enumerate(distances):
            if i != 0:
                groups[i] = groups[i-1]
                if distances[i]-distances[i-1] > h:
                    groups[i] += 1
                    group = groups[i]
        averaged_distances = nu.array([nu.sum(distances[nu.where(groups == i)])/len(nu.where(groups == i)[0]) for i in range(group+1)])
        epsilon = nu.zeros((N,N), dtype=int)
        for (i,j), g in zip(index_list, groups):
            epsilon[i,j] = g
        mask = nu.zeros(nu.shape(epsilon), dtype=int)
        for i in range(N):
            for j in range(N):
                if i != j and averaged_distances[epsilon[i,j]] < cut_radius:
                    mask[i,j] = 1
        averaged_distances = averaged_distances[nu.where(averaged_distances < cut_radius)]
        return epsilon, averaged_distances, mask


    def _find_forces(self, function, v_rep_points):
        """
        Try to minimize the residual forces by finding matrix
        elements epsilon_ij = V'_rep(|r_ij|).
        """
        from scipy.optimize import fmin
        last_res_forces = 0
        N = len(v_rep_points)
        while True:
            ret = fmin(function, v_rep_points, full_output=1, disp=0)
            v_rep_points = ret[0]
            forces = ret[1]
            if ret[4] == 0 and nu.abs(forces-last_res_forces) < 0.001:
                return v_rep_points, forces, True
            last_res_forces = ret[1]
        return v_rep_points, last_res_forces, False


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


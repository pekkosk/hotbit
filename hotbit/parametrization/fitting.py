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

class RepulsiveFitting:

    def __init__(self,symbol1,symbol2,r_dimer=None,calc=None,maxiter=None,errfile=None): 
        """
        Fit the short-range repulsive potential.
        
        Parameters:
        -----------
        rep:
        r_cut:
        r_dimer:
        order:
        calc:
        """
        #raise NotImplementedError('Not yet consistently implemented; however, you may comment this exception to play around.')
        self.elm1=Element(symbol1)
        self.elm2=Element(symbol2)
        self.sym1=symbol1
        self.sym2=symbol2
        self.r_dimer=r_dimer
        
        self.deriv=[]
        self.comments=''
        self.scale=1.025                    # scaling factor for scalable systems
        self.structures = []
        self.calc=calc
        self.v=None
        self.r_cut=None                    # the cutoff radius, set in fitting
        self.maxiter=maxiter

        self.colors = ['blue','cyan','red','pink','yellow','orange','#8DEE1E','magenta','green','white','black']
        self.color_index = 0

        self.err = None
        self._set_err_out(errfile)
        print '\n\n\n\nCreating a repulsion curve between elements %s and %s' % (self.sym1, self.sym2)
        print '  r_dimer      =',r_dimer
        print '  1.5 x r_dimer=',1.5*r_dimer


    def __call__(self,r,der=0):
        """ Return V_rep(r) or V_rep'(r) """
        return self.v(r,der=der)


    def _set_err_out(self, err):
        import sys
        if self.err not in [None, sys.stdout]:
            if hasattr(self.err, 'close'):
                try:
                    self.err.close()
                except:
                    pass
                self.err = None
        if type(err) == str:
            self.err = open(err, 'w')
        elif hasattr(err, 'write'):
            self.err = out
        else:
            self.err = sys.stdout


    def plot(self, filename=None):
        """ Plot vrep and derivative together with fit info. 
        
        parameters:
        ===========
        filename:     graphics file name
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
        pl.ylim(ymin=0,ymax=self(0))
        pl.xlim(xmin=0, xmax=self.r_cut)

        # Vrep'
        pl.subplot(1,2,2)
        pl.ylabel(r'$dV_{rep}(r)/dr$ (eV/$\AA$)')
        pl.xlabel(r'$r$ ($\AA$)')
        pl.plot(r,vp,label=r'$dV_{rep}(r)/dr$')
        for s in self.deriv:
            pl.scatter( [s[0]],[s[1]],s=100*s[2],c=s[3],label=s[4])
        pl.axvline(x=self.r_cut,c='r',ls=':')
        pl.axvline(x=self.r_dimer,c='r',ls=':')
        xmin = 0.8*self.r_dimer
        xmax = 1.2*self.r_cut
        ymin = self(xmin, der=1)
        ymax = 0
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
        """ Append some comment for par-file. """
        if s==None: return
        if s in [None, '']:
            return
        add='|'
        if len(self.comments)==0:
            add=''
        self.comments+=add+s


    def _scale_positions(self,x,cell=False):
        """ Scale the whole system by x; also the unit cell. """
        self.set_cell(self.get_cell()*x,fix=False)
        self.reduce_atoms_into_cell()


    def _solve_ground_state(self, atoms, charge=0, calc=None):
        """
        Solves the ground state of given atoms and returns the
        solved calculator.
                
        """
        if calc == None:
            c = copy(self.calc)
            c.set('txt', self.calc.txt)
        else:
            c = copy(calc)
            c.set('txt', calc.txt)
        c.set("charge",charge)
        atoms.set_calculator(c)
        atoms.get_potential_energy()
        c.set('txt',"-")
        return c


    def fit(self, r_cut, s=None): #, filename=None):
        """
        Fit smoothing spline into points {r, V_rep'(r)}.
        The weights are treated as inverse of the stardard deviation.
        """
        from scipy.interpolate import splrep, splev
        
        # fitting parameters; order of spline is cubic at maximum
        k = min(len(self.deriv)-1,3)
        self.r_cut = r_cut
        self.s = s
        self.k = k
        
        
        x = nu.array([self.deriv[i][0] for i in range(len(self.deriv))])
        y = nu.array([self.deriv[i][1] for i in range(len(self.deriv))])
        w = nu.array([self.deriv[i][2] for i in range(len(self.deriv))])
        # sort values so that x is in ascending order
        indices = x.argsort()
        x = x[indices]
        y = y[indices]
        w = w[indices]
        x, y, w = self._unite_closeby_points(x,y,w)
        # use only points that are closer than r_cut
        indices = nu.where(x < r_cut)
        x = list(x[indices])
        y = list(y[indices])
        w = list(w[indices])
        # force the spline curve to go to zero at x=r_cut
        x.append(r_cut)
        y.append(0)
        w.append(1e3*max(w))
        if s == None:
            # from documentation of splrep in scipy.interpolate.fitpack
            s = len(x) - nu.sqrt(2*len(x))
        print ""
        print "  Fitting smoothing spline with parameters"
        print "  k=%i, s=%0.4f, r_cut=%0.4f\n" %(k, s, r_cut)
        tck = splrep(x, y, w, s=s, k=k)

        def dv_rep(r):
            return splev(r, tck)

        v_rep = self._integrate_vrep(dv_rep, r_cut)

        def potential(r, der=0):
            if der == 0:
                return v_rep(r)
            elif der == 1:
                return dv_rep(r)
            else:
                raise NotImplementedError("Only 0th and 1st derivatives")
        self.v = potential
        

    def _unite_closeby_points(self, x, y, w):
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


    def write_par(self, filename=None):
        """
        Write the full par-file to file.
        
        parameters:
        ===========
        filename:     output file
        """
        from time import asctime
        import shutil
        r_cut, s, k = self.r_cut, self.s, self.k
        par_files = self.calc.tables
        par_file = None
        e12 = self.sym1 + self.sym2
        e21 = self.sym2 + self.sym1
        if e12 in par_files:
            par_file = par_files[e12]
        if e21 in par_files:
            par_file = par_files[e21]
        if par_file != None:
            shutil.copy(par_file, filename)
            f = open(filename, 'a')
            # add comments
            print >> f, "repulsion_comment="
            print >> f, "%s\nparameters r_cut = %0.4f Ang, s = %0.4f, k = %3i" % (asctime(),r_cut, s, k)
            if len(self.structures)>1:
                print >> f, "The systems used to produce this fit:"
                for data in self.structures:
                    print >> f, "%20s %3s" % (data['filename'], data['charge'])
            if len(self.comments) > 0:
                print >> f, self.comments            
            print >> f, '\n\nrepulsion='
            for r in nu.linspace(0.1, r_cut, 100):
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
        # SplineFunction wants the x-values in ascendind order
        return SplineFunction(r_g[::-1], v_rep[::-1])


    #
    #       Fitting methods
    #             
    def append_point(self,weight,r,dvrep,comment=None,label='point',color='g'):
        """ Add point to vrep'-fitting.
        
        parameters:
        ===========
        weight:    relative weight (inverse of deviation)
        r:         radius, in Angstroms
        vder:      value of the derivative at r
        color:     pylab color for plotting
        info:      label in plotting
        """
        self.deriv.append([r,dvrep,weight,color,label])
        self.add_comment(comment)


    def append_dimer(self,weight,r_dimer=None,comment=None,label='dimer',color='r'):
        """ Use dimer bond length in fitting. 
        
        parameters:
        ===========
        weight:    relative weight (inverse of deviation)
        r_dimer:   dimer bond length in Angstroms. By default the initialization value.        
        """
        d=r_dimer
        if d==None: d=self.r_dimer
        dimer=Atoms(symbols=[self.sym1,self.sym2],\
                    positions=[(0,0,0),(d,0,0)],\
                    pbc=False,cell=[100,100,100])
        print dimer
        self.append_scalable_system(weight,dimer,0.0,comment=comment,label=label,color=color)


    def append_scalable_system(self,weight,system,charge,comment=None,label=None,color='m'):
        """ Use scalable equilibrium (DFT) system in repulsion fitting. """
        #raise NotImplementedError('Not implemented correctly')
        if type(system)==type(''):
            from ase import read
            atoms=Atoms(read(system))
            atoms.center(vacuum=10)
            if name == None:
                name=system+'_scalable'
        else:
            atoms=system
            if label == None:
                label=atoms.get_name()+'_scalable'
        #forces = self.calc.get_forces(atoms)
        #for vec_f in forces:
        #    if nu.linalg.norm(vec_f) > 0.05:
        #        raise Exception("System is not in an equilibrium!")
        x=self.scale
        r=atoms.mean_bond_length()
        bonds=atoms.number_of_bonds()
        e1=self._solve_ground_state(atoms,charge).get_potential_energy(atoms)
        atoms.scale_positions(x)
        e2=self._solve_ground_state(atoms,charge).get_potential_energy(atoms)

        dEdr=(e2-e1)/(x*r-r)
        self.append_point(weight,r,-dEdr/bonds,comment,label,color)
        if comment!=None:
            self.add_comment(comment)


    def append_energy_curve(self, weight, trajectory, charge=0, comment=None, cutoff=3, toldist=1e-5, label='energy curve', color='y'):
        """
        Calculates the V'rep(r) from a given DFT ase-trajectory for elements
        A and B:

                    E'_DFT(R) - E'_BS(R)
        V'rep(R) =  ------------------ ,
                            N

        where R is the distance between the elements in the system and
              N is the number of different A-B pairs that are taken
              into account.

        parameters:
        ===========        
        weight:              Fitting weight (inverse of deviation)
        trajectory:          ase trajectory for energy curve
        charge:              total charge of the system
        comment:             A comment to be added to the .par file
        cutoff:              Cutoff for element interaction, in WHAT UNITS?
        toldist:             Acceptable tolerance in distance
                             (such that R and R+toldist still have the 'same' bond length)        
        label:               Name that is used in the fitting plot.
        color:               Color that is used in the fitting plot.
        """
        import scipy
        from ase.io.trajectory import PickleTrajectory
        from ase import Atoms
        from box.interpolation import SplineFunction

        sigma = 1.0/weight
        calc = self.calc
        if color == None:
            color = self._get_color()
        print ""
        print "  Appending energy curve data from %s..." % trajectory
        traj = PickleTrajectory(trajectory)
        R, E_dft, N = self._process_trajectory(traj, self.sym1, self.sym2, cutoff, toldist)
        E_bs = nu.zeros(len(E_dft))
        usable_frames = []
        for i in range(len(traj)):
            atoms=copy(traj[i])
            calc_new = self._solve_ground_state(atoms, charge, calc)
            if calc_new == None:
                print "    *** Error: No data from frame %i ***" % i
                print >> self.err, "No data from %s frame %i" % (trajectory, i)
            else:
                E_bs[i] = calc_new.get_potential_energy(atoms)
                usable_frames.append(i)
                calc_new.timer.summary()
                calc_new.set_text("-")
                del(calc_new)
                print "    Collected data from frame %i" % i
        traj.close()
        # use only frames where DFTB-calculation converged
        R = R[usable_frames]
        E_dft = E_dft[usable_frames]
        E_bs = E_bs[usable_frames]
        sigma_coefficient = nu.sqrt(len(R))
        # sort the radii and corresponding energies to ascending order
        indices = R.argsort()
        R = R[indices]
        E_dft = E_dft[indices]
        E_bs = E_bs[indices]
        vrep = SplineFunction(R, (E_dft - E_bs)/N)

        for i, r in enumerate(R):
            if i > 0:
                label='_nolegend_'
            self.append_point(1./(sigma_coefficient*sigma),r, vrep(r,der=1), comment, label, color)


    def append_homogeneous_structure(self, weight, filename, charge=0, comment=None, cutoff=3, toldist=0.005, traj_indices=None,label='homogeneous structure',color=None):
        """
        For a given structure, calculate points {r, V'_rep(r)} so that
        the residual forces are minimized (F_i = \sum_j(dV(|r_ij|)/dR)).
        If only coordinates are given, the structure must be an equilibrium
        structure. If also forces from DFT calculation are given
        (ase.traj), any homogeneous structure may be applied, assuming
        the minimization converges.
        
        parameters:
        ===========
        weight:        relative weight
        filename:      input structure
        charge:        structure's electric charge
        comment:       fitting comment
        cutoff:        interaction cutoff
        toldist:       acceptable tolerance in bond lengths (still considered equal)
        traj_indices:  A list of indices of the frames that is read from
                       the trajectory.
        label:         plotting label
        color:         plotting color
        """
        sigma = 1.0/weight
        points = []
        structures, traj_indices = self.import_structures(filename, traj_indices)
        print ""
        print "Appending homogeneous structure from %s..." % filename
        for ind, (fr, structure) in enumerate(zip(traj_indices, structures)):
            print "  frame %i" % fr
            epsilon, distances, mask=self._get_matrices(structure, cutoff, toldist)
            if len(distances) == 0:
                raise RuntimeError("There are no bonds under given cut radius in frame %i." & fr)
            r_hat = self._construct_rhat_matrix(structure)

            # if the given structure contains forces, use them, otherwise
            # assume that the structure is optimized
            N = len(structure)
            try:
                forces_DFT = structure.get_forces()
                print "    Found forces from trajectory %s at frame %i" % (filename, fr)
            except:
                forces_DFT = nu.zeros((N,3))
                print "    No forces in trajectory %s at frame %i" % (filename, fr)
            calc_new = self._solve_ground_state(structure, charge=charge)
            if calc_new == None:
                print >> self.err, "    No data from %s frame %i" % (label, fr)
            else:
                forces = calc_new.get_forces(structure)
                forces_res = forces_DFT - forces

                # use one less point in an array that is given for
                # the minimizer to reduce the number of degrees of freedom
                v_rep_points = nu.zeros(len(distances)-1)

                # the function that is minimized, returns the sum of
                # the norm of forces acting on each atoms
                def residual_forces(v_rep_points):
                    res = 0.
                    for i in range(N):
                        indices = nu.array(mask[i,:]*epsilon[i,:])
                        indices -= 1
                        ppoints  = v_rep_points[indices]
                        f = forces_res[i].copy()
                        f -= nu.dot((mask[i,:]*ppoints),r_hat[i,:])
                        res += nu.dot(f,f)
                    return res

                print "    Found %i different bond lengths." % len(v_rep_points)
                v_rep_points, last_res_forces, minimized = self._find_forces(residual_forces, v_rep_points)
                print "    The sum of the squared norms of the net forces: %0.4f (eV/ang)**2" % (last_res_forces)
                # finally add the missing component
                v_rep_points = list(v_rep_points)
                v_rep_points.insert(0,0)
                if minimized:
                    for i in range(0,N-1):
                        for j in range(i+1,N):
                            if mask[i,j] > 0:
                                # could add a degeneracy factor since there
                                # may be different number of different bonds
                                points.append([distances[epsilon[i,j]], v_rep_points[epsilon[i,j]], last_res_forces])
                else:
                    print "    The minimization of forces did not converge!"
        if len(points) > 0:
            points = nu.array(points)
            sigmas = points[:,2] + nu.min(points[:,2])*0.001 # add small value to prevent 1/sigma to go infinity
            inv_sigmas = 1./sigmas
            norm_factor = nu.sqrt(1./sigma**2 / nu.dot(inv_sigmas, inv_sigmas))
            points[:,2] /= norm_factor
            if color == None:
                color = self._get_color()
            for data in points:
                self.append_point(1/data[2], data[0], data[1], comment, label, color)
                label = '_nolegend_'
            self.add_comment(comment)

    def _process_trajectory(self, traj, elA, elB, separating_distance, h):
        """
        Check each frame in trajectory, make sure that the trajectory
        is suitable for repulsion fitting for the elements A and B.
        Finally returns the bond lengths of elements A and B and
        the DFT energy in each image.
        """

        E_dft = nu.zeros(len(traj))
        R = nu.zeros(len(traj))
        N = nu.zeros(len(traj))
        self._assert_fixed_bond_lengths_except(traj, elA, elB, separating_distance, h)
        for i, image in enumerate(traj):
            atoms = copy(Atoms(image))
            r, n = self._get_distance_of_elements(elA, elB, atoms, separating_distance, h)
            E_dft[i] = image.get_total_energy()
            R[i] = r
            N[i] = n
        assert nu.all(N==N[0])
        return R, E_dft, N[0]


    def _assert_fixed_bond_lengths_except(self, t, elA, elB, separating_distance, h):
        """
        Makes sure that all the element pairs expect pairs A-B are the same 
        or larger than the defined limit in all configurations.
        """
        fixed_pairs = []
        fixed_lengths = []
        long_pairs = []
        atoms = t[0]
        for i in range(len(atoms)-1):
            for j in range(i+1,len(atoms)):
                a = atoms[i]
                b = atoms[j]
                dL = nu.linalg.norm(a.position-b.position)
                if not (( a.symbol == elA and b.symbol == elB ) or \
                       ( a.symbol == elB and b.symbol == elA )):
                    if dL > separating_distance:
                        long_pairs.append([i,j])
                    else:
                        fixed_pairs.append([i,j])
                        fixed_lengths.append(dL)
        for k in range(1, len(t)):
            atoms = t[k]
            for i in range(len(atoms)-1):
                for j in range(i+1,len(atoms)):
                    a = atoms[i]
                    b = atoms[j]
                    dL = nu.linalg.norm(a.position-b.position)
                    if [i,j] in fixed_pairs:
                        index = fixed_pairs.index([i,j])
                        if nu.abs(dL - fixed_lengths[index]) > h:
                            raise AssertionError("Fixed bond length vary too much in the trajectory: atoms %i and %i." % (i, j))
                    if [i,j] in long_pairs:
                        if dL < separating_distance:
                            raise AssertionError("Long bond goes below the separating limit: atoms %i and %i." % (i, j))


    def _get_distance_of_elements(self, elA, elB, positions, separating_distance, h):
        """
        Calculates the distances of element pairs A-B that are closer to
        each other than the defined limit and makes sure that the distances
        are equal.
        """
        R = []
        for i in range(len(positions)-1):
            for j in range(i+1, len(positions)):
                a = positions[i]
                b = positions[j]
                if ( a.symbol == elA and b.symbol == elB ) or \
                   ( a.symbol == elB and b.symbol == elA ):
                    R.append(nu.linalg.norm(a.position - b.position))
        R.sort()
        R = nu.array(R)
        R_min = R[0]
        N = nu.sum(nu.where(nu.abs(R - R_min) < h, 1, 0))
        for r in R[N:]:
            if r < separating_distance:
                raise AssertionError("Element pairs are too close to each other.")
        return nu.average(R[:N]), N


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


    def _get_matrices(self, structure, cut_radius, h):
        """
        Construct epsilon matrix that maps the indices (i,j) to a
        single list of distances. If there are many bonds with
        almost same lengths, treat these bonds as there was only
        one of them in order to reduce the degree of freedom
        in the minimization. If the bond is longer that given
        cut radius, that bond is ignored.
        """
        N = len(structure)
        distances = []
        index_list = []
        for i in range(N):
            for j in range(N):
                vec = structure.positions[j] - structure.positions[i]
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


    def append_homogeneous_bulk(self, weight, trajectory, coordination, comment=None, cutoff=3.0, toldist=1e-5, label='bulk', color=None):
        #raise Exception("Not well tested, comment this line to use this method")
        """ Get repulsion fitting from homogeneous bulk calculation.
            The coordination number is the no. of nearest-neighbours
            per atom in the bulk (e.g 4 in diamond lattice).
            
            
        parameters:
        ===========
        weight:        fitting weight
        trajectory:
        coordination: 
        comment:       comment fitting
        cutoff:        interaction cutoff
        toldist:       tolerance for distances (still considered the same)
        label:         plotting label
        color:         plotting color
        """
        print ""
        print "Appending homogeneous bulk from %s..." % trajectory
        sigma=1.0/weight
        if type(trajectory) == str:
            traj = PickleTrajectory(trajectory)
        for image in traj:
            R, E_dft, N = self._process_trajectory(traj, self.sym1, self.sym2, cutoff, toldist)
        N_a = len(traj[0])
        N = N_a * coordination / 2.

        E_bs = nu.zeros(len(E_dft))
        usable_frames = []
        for i in range(len(traj)):
            atoms=copy(traj[i])
            calc_new = self._solve_ground_state(atoms, 0.0, self.calc)
            if calc_new == None:
                print "    *** Error: No data from frame %i ***" % i
                print >> self.err, "No data from %s frame %i" % (dft_traj, i)
            else:
                E_bs[i] = calc_new.get_potential_energy(atoms)
                usable_frames.append(i)
                calc_new.timer.summary()
                calc_new.set_text("-")
                del(calc_new)
                print "    Collected data from frame %i" % i
        traj.close()
        # use only frames where DFTB-calculation converged
        R = R[usable_frames]
        E_dft = E_dft[usable_frames]
        E_bs = E_bs[usable_frames]
        sigma_coefficient = nu.sqrt(len(R))
        # sort the radii and corresponding energies to ascending order
        indices = R.argsort()
        R = R[indices]
        E_dft = E_dft[indices]
        E_bs = E_bs[indices]
        from box.interpolation import SplineFunction
        vrep = SplineFunction(R, (E_dft - E_bs)/N)
        if color == None:
            color = self._get_color()

        for i, r in enumerate(R):
            if i > 0:
                label='_nolegend_'
            self.append_point(1./(sigma_coefficient*sigma), r, vrep(r, der=1), comment, label, color)


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


    def import_structures(self, filename, traj_indices=None):
        structures = []
        if ".traj" in filename:
            traj = PickleTrajectory(filename)
            if traj_indices == None:
                traj_indices = range(len(traj))
            for i in traj_indices:
                image = traj[i]
                structures.append(image)
        elif ".xyz" in filename:
            structure = read(filename)
            structure.center(vacuum=6)
            structures.append(structure)
            traj_indices = [0]
        else:
            raise Exception("Unknown file format: %s" % structure)
        return structures, traj_indices


    def _get_color(self):
        color = self.colors[self.color_index]
        self.color_index += 1
        if self.color_index == len(self.colors):
            self.color_index = 0
        return color


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
                print ex
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


from hotbit import Element
from hotbit import Calculator
import numpy as nu
from ase.units import Bohr,Hartree
from box import mix
from box import Atoms
from box.interpolation import Function
from ase import read, PickleTrajectory
import pylab as pl


class RepulsiveFitting:
    
    def __init__(self,rep,r_cut=None,r_dimer=None,order=8,calc=None,maxiter=None):
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
        self.elm1=Element(rep[0])
        self.elm2=Element(rep[1])
        self.sym1=self.elm1.get_symbol()
        self.sym2=self.elm2.get_symbol()
        self.r_dimer=r_dimer
        self.r_cut=r_cut                    # the cutoff radius
        self.r_small=r_dimer*0.5            # use ZBL repulsion for r<r_small
        self.order=order                    # order of Vrep polynomial 
        self.deriv=[]
        self.comments=''
        self.scale=1.025                    # scaling factor for scalable systems
        self.calc=calc
        self.v=None
        self.maxiter=maxiter

        self.colors = ['blue','cyan','red','pink','yellow','orange','#8DEE1E','magenta','green','white','black']
        self.color_index = 0

        print 'r_dimer      =',r_dimer
        print '1.5 x r_dimer=',1.5*r_dimer


    def __call__(self,r,der=0):
        """ Return V_rep(r) or V_rep'(r) """
        if self.v is None:
            return self.parametric_potential(r,self.param,der=der)
        else:
            return self.v(r,der=der)


    def write_to_par(self,txt=None,points=100,append=True):
        """ Append Vrep into par-file, including comments. 
        
        POISTETAAN. Liittyy fortran-versioon.
        """
        if txt is None:
            txt='%s_%s.par' %(self.elm1,self.elm2)
        mode=('w','a')[append]
        o=open(txt,mode)
        print>>o, 'fitting='
        print>>o, self.comments
        print>>o, '\n\nrepulsion='
        for r in nu.linspace(0.1,self.r_cut,points):
            print>>o, r/Bohr, self(r)/Hartree
        o.close()        
        
        
    def plot(self):
        """ Plot vrep and derivative together with fit info. """
        r=nu.linspace(self.r_small,self.r_cut)
        v=[self(x,der=0) for x in r]
        vp=[self(x,der=1) for x in r]
        rmin=0.9*min([d[0] for d in self.deriv])

        # Vrep
        pl.subplots_adjust(wspace=0.25)
        pl.subplot(1,2,1)
        pl.ylabel('$V_{rep}(r) (eV)$')
        pl.xlabel('$r (\AA)$')
        pl.axvline(x=self.r_cut,c='r',ls=':')
        pl.plot(r,v)
        pl.ylim(ymin=0,ymax=self(rmin))
        pl.xlim(xmin=rmin)

        # Vrep'
        pl.subplot(1,2,2)
        pl.ylabel('$dV_{rep}(r)/dr (eV/\AA)$')
        pl.xlabel('$r (\AA)$')
        pl.plot(r,vp,label='$dV_{rep}(r)/dr$')
        for s in self.deriv:
            pl.scatter( [s[0]],[s[1]],s=100*s[2],c=s[3],label=s[4])
        pl.axvline(x=self.r_cut,c='r',ls=':')
        pl.axvline(x=self.r_dimer,c='r',ls=':')
        xmin = 0.8*self.r_dimer
        xmax = 1.2*self.r_cut
        ymin = self(0.8*self.r_dimer, der=1)
        ymax = abs(0.2*ymin)
        pl.text(self.r_dimer, ymin, 'r_dimer')
        pl.text(self.r_cut, ymin, 'r_cut')
        pl.xlim(xmin=xmin, xmax=xmax)
        pl.ylim(ymin=ymin, ymax=ymax)
        pl.legend()
        pl.show()


    def add_fitting_comment(self,s):
        """ Append some comment for par-file. """
        add='|'
        if len(self.comments)==0:
            add=''
        self.comments+=add+s    
        
            
    def ZBL_derivative(Z1,Z2,r):
        """ 
        Return the derivative of the Ziegler,Biersack,Littmar potential
        for elements Z1 and Z2 at distance r.
        """
        def phi(x,der=0):
            if der==0:
                return 0.1818*exp(-3.2*x) + 0.5099*exp(-0.9423*x) + 0.2802*exp(-0.4029*x) \
                    +0.02817*exp(-0.2016*x)
            elif der==1:    
                return -0.58176*exp(-3.2*x) - 0.48047877*exp(-0.9423*x) - \
                    0.11289258*exp(-0.4029*x) - 0.005679072*exp(-0.2016*x)
        au=0.8854*0.5292/( Z1**0.23+Z2**0.23 )
        return -Z1*Z2/r**2*phi(r/au) + Z1*Z2/r*phi(r/au,der=1)/au
            
            
    def parametric_potential(self,r,d,der=0):
        """ 
        Return Vrep(r) with given parameter representation. 
        
        Vrep(r) =sum_i[0,order] d_i (r_c-r)**i
        Vrep'(r)=sum_i[1,order] d_i (r_c-r)**(i-1)*i*(-1)
        """
        from copy import copy
        re=0.0
        d2=copy(d)
        d2[0:2]=0.0 # first two terms have to be zero
        if r>self.r_cut or r<0:
            return 0.0
        else:
            if der==0:
                return sum( [d2[i]*(self.r_cut-r)**i for i in range(self.order)] )
            elif der==1:
                return sum( [d2[i]*(self.r_cut-r)**(i-1)*i*(-1) for i in range(1,self.order)] )


    def scale_positions(self,x,cell=False):
        """ Scale the whole system by x; also the unit cell. """
        self.set_cell(self.get_cell()*x,fix=False)
        self.reduce_atoms_into_cell()
        
        
    def solve_ground_state(self, atoms, charge=None, calc=None):
        """
        vahan kuten get_energy
                
        """
        from copy import copy
        if calc == None:
            c = copy(self.calc)
        else:
            c = copy(calc)
        if charge != None:
            c.set("charge",charge)
        atoms.set_calculator(c)
        try:
            # FIXME make these output to file also
            atoms.get_potential_energy()
        except Exception:
            del(c)
            c = None
        return c


    def get_energy(self,atoms,charge,forces=False):
        """
        Calculate energy (or forces) for given structure with given charge
        and given parameters (elements & tables as well).
        
        MUUTETAAN.(poistetaan?)
        """
        from copy import deepcopy, copy
        if type(atoms)==type(''):
            atms=read(file)
        else:
            atms=atoms
        calc=copy(self.calc)
        calc.set(charge=charge)
        atms.set_calculator(calc)
        if forces:
            res=atms.get_forces()
        else:
            res=atms.get_potential_energy()
        #calc.finalize()
        return res



    def fitting_function(self,d):
        """ Minimize this function in the fitting. """
        chi2=0.0
        for point in self.deriv:
            r=point[0]
            vp=point[1]
            w=point[2]
            chi2+=(vp-self.parametric_potential(r,d,der=1))**2*w
        return chi2
        

    def fit(self, r_cut):
        """ Fit V_rep(r) into points {r,V_rep'(r)}. """
        from scipy.optimize import fmin
        self.r_cut = r_cut
        self.param=nu.zeros(self.order,float)
        self.param[2]=10                    # initial guess...
        print "Fitting with r_cut=%0.6f..." % r_cut
        self.param = fmin(self.fitting_function,self.param,maxiter=self.maxiter, maxfun=self.maxiter)


    def fit_smoothing_spline(self, r_cut, s=None, filename=None):
        """
        Fit smoothing spline into points {r, V_rep'(r)}.
        The weights are treated as inverse of the stardard deviation.
        """
        from scipy.interpolate import splrep, splev
        self.r_cut = r_cut
        k = 3
        x = nu.array([self.deriv[i][0] for i in range(len(self.deriv))])
        y = nu.array([self.deriv[i][1] for i in range(len(self.deriv))])
        w = nu.array([self.deriv[i][2] for i in range(len(self.deriv))])
        # sort values so that x is in ascending order
        indices = x.argsort()
        x = x[indices]
        y = y[indices]
        w = w[indices]
        x, y, w = self.average_too_similar_values(x,y,w)
        # use only points that are closer than r_cut
        indices = nu.where(x<r_cut)
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
        print "Fitting smoothing spline with parameters"
        print "k=%i, s=%i, r_cut=%0.4f\n" %(k, s, r_cut)
        tck = splrep(x, y, w, s=s, k=k)

        def dv_rep(r):
            return splev(r, tck)

        v_rep = self.integrate_vrep(dv_rep, r_cut)

        def potential(r, der=0):
            if der == 0:
                return v_rep(r)
            elif der == 1:
                return dv_rep(r)
            else:
                raise NotImplementedError("Only 0th and 1st derivatives")
        self.v = potential
        if filename != None:
            self.smoothing_spline_fit_to_file(filename, r_cut, s, k)


    def average_too_similar_values(self, x, y, w):
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


    def smoothing_spline_fit_to_file(self, filename, r_cut, s, k):
        """
        Write the full par-file to filename.
        """
        from time import asctime
        import shutil
        par_files = self.calc.table_files
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
            print >> f, "repulsion_fitting_comment="
            print >> f, "Repulsive potential generated by fitting function 'fit_smoothing_spline'"
            print >> f, "parameters r_cut = %0.4f, s = %0.4f, k = %3i" % (r_cut, s, k)
            print >> f, asctime() + "\n\n"
            print >> f, 'fitting='
            print >> f, self.comments
            print >> f, '\n\nrepulsion='
            for r in nu.linspace(0.1, r_cut, 100):
                print >> f, r/Bohr, self(r)/Hartree
            f.close()


    def integrate_vrep(self, dv_rep, r_cut, N=100):
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


    def repulsion_forces(self,atoms,vrep):
        """ 
        Return the repulsive forces for atoms using given
        vrep(r) function for present element pairs.
        
        POISTETAAN: kayta calc:n 
        calc.rep.get_repulsive_forces()
        """
        raise NotImplementedError('No idea if this works, check!')
        
        n=len(atoms)
        forces=nu.empty((n,3))
        els=atoms.get_chemical_symbols()
        
        for i in range(n):
            for j in range(i,n):
                if (els[i],els[j]) is (self.elm1,self.elm2) or \
                   (els[j],els[i]) is (self.elm1,self.elm2):
                    rij=atoms.vector(i,j)
                    r=nu.linalg.norm(rij)
                    force=vrep(r,der=1)*rij/r
                    forces[i,:]=force
                    forces[j,:]=-force
        return forces
                    
    #
    #       Fitting methods
    #             
    def append_point(self,data,comment=None):
        """ Add point to vrep'-fitting: data=[r,v',w,color,info] """
        self.deriv.append(data)
        if comment is not None:
            self.add_fitting_comment(comment)
           
           
    def append_dimer(self,weight):
        """ Use dimer bond length in fitting. """
        dimer=Atoms(symbols=[self.sym1,self.sym2],\
                    positions=[(0,0,0),(self.r_dimer,0,0)],\
                    pbc=False,cell=[100,100,100])
        self.append_scalable_system(dimer,0.0,weight,name='dimer bond length',comment='dimer at %.4f %s' %(self.r_dimer,chr(197)))


    def append_scalable_system(self,system,charge,weight,color='#71FB00',name=None,comment=None):
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
            if name == None:
                name=atoms.get_name()+'_scalable'
        #forces = self.calc.get_forces(atoms)
        #for vec_f in forces:
        #    if nu.linalg.norm(vec_f) > 0.05:
        #        raise Exception("System is not in an equilibrium!")
        x=self.scale
        r=atoms.mean_bond_length()
        bonds=atoms.number_of_bonds()
        e1=self.solve_ground_state(atoms,charge).get_potential_energy(atoms)
        atoms.scale_positions(x)
        e2=self.solve_ground_state(atoms,charge).get_potential_energy(atoms)

        dEdr=(e2-e1)/(x*r-r)
        self.append_point([r,-dEdr/bonds,weight,color,name])
        if comment is None:
            comment='scalable %s' %name
        self.add_fitting_comment(comment)


    def append_energy_curve(self, dft_traj, **kwargs):
        """
        Calculates the V'rep(r) from a given DFT ase-trajectory for elements
        A and B:

                    E'_DFT(R) - E'_BS(R)
        V'rep(R) =  ------------------ ,
                            N

        where R is the distance between the elements in the system and
              N is the number of different A-B pairs that are taken
              into account.

        acceptable keyword arguments:
        color:               Color that is used in the fitting plot.
        label:               Name that is used in the fitting plot.
        sigma:               The "uncertainty" of the points calculated
                             from this system (default=1).
        separating_distance: If distance between two elements is larger
                             than this, it is assumed that their
                             interaction can be neglected.
        h:                   The deviation that is allowed in the bond
                             length between the element pairs that should
                             have exactly the same bond length.
        comment:             A comment to be added to the .par file
        """
        import scipy
        from ase.io.trajectory import PickleTrajectory
        from ase import Atoms
        from copy import copy
        from box.interpolation import SplineFunction

        if not 'separating_distance' in kwargs:
            kwargs['separating_distance'] = 3.0
        if not 'h' in kwargs:
            kwargs['h'] = 1e-5
        if not 'sigma' in kwargs:
            kwargs['sigma'] = 1.0
        if not 'charge' in kwargs:
            kwargs['charge'] = None
        traj = PickleTrajectory(dft_traj)
        R, E_dft, N = self.process_trajectory(traj, self.sym1, self.sym2, **kwargs)
        E_bs = nu.zeros(len(E_dft))
        usable_frames = []
        for i in range(len(traj)):
            atoms=copy(traj[i])
            calc = self.solve_ground_state(atoms, kwargs['charge'])
            if calc != None:
                E_bs[i] = calc.get_potential_energy(atoms)
                del(calc)
                usable_frames.append(i)
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

        if not 'color' in kwargs:
            color = self.get_color()
        else:
            color = kwargs['color']
        if not 'label' in kwargs:
            label = ''
        else:
            label = kwargs['label']
        if not 'comment' in kwargs:
            comment = "Point from energy curve fitting"
        else:
            comment = kwargs['comment']
        for i, r in enumerate(R):
            if i > 0:
                label='_nolegend_'
            self.append_point([r, vrep(r,der=1), 1./(sigma_coefficient*kwargs['sigma']), color, label], comment=comment)


    def process_trajectory(self, traj, elA, elB, **kwargs):
        """
        Check each frame in trajectory, make sure that the trajectory
        is suitable for repulsion fitting for the elements A and B.
        Finally returns the bond lengths of elements A and B and
        the DFT energy in each image.
        """
        from copy import copy

        E_dft = nu.zeros(len(traj))
        R = nu.zeros(len(traj))
        N = nu.zeros(len(traj))
        self.assert_fixed_bond_lengths_except(traj, elA, elB, **kwargs)
        for i, image in enumerate(traj):
            atoms = copy(Atoms(image))
            r, n = self.get_distance_of_elements(elA, elB, atoms, **kwargs)
            E_dft[i] = image.get_total_energy()
            R[i] = r
            N[i] = n
        assert nu.all(N==N[0])
        return R, E_dft, N[0]


    def assert_fixed_bond_lengths_except(self, t, elA, elB, **kwargs):
        """
        Makes sure that all the element pairs expect pairs A-B are the same 
        or larger than the defined limit in all configurations.
        """
        separating_distance = kwargs['separating_distance']
        h = kwargs['h']

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


    def get_distance_of_elements(self, elA, elB, positions, **kwargs):
        """
        Calculates the distances of element pairs A-B that are closer to
        each other than the defined limit and makes sure that the distances
        are equal.
        """
        separating_distance = kwargs['separating_distance']
        h = kwargs['h']

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


    def append_homogeneous_structure(self, filename, charge=0, color=None, weight=1.0, label='homogeneous structure', comment='', cut_radius=3, traj_indices=None, h=0.005, calc=None):
        """
        For a given structure, calculate points {r, V'_rep(r)} so that
        the residual forces are minimized (F_i = \sum_j(dV(|r_ij|)/dR)).
        If only coordinates are given, the structure must be an equilibrium
        structure. If also forces from DFT calculation are given
        (ase.traj), any homogeneous structure may be applied, assuming
        the minimization converges.

        maxiter:      how many fmin-runs is performed to find minimum point
        cut_radius:   the largest distance of elements that is taken
                      into account.
        traj_indices: list of indices for images in trajectory
        h:            the variation in the bond lengths that are still
                      considered to be equal.
        """
        #raise NotImplementedError("Not tested adequately")
        points = []
        structures = self.import_structures(filename, traj_indices)
        for structure in structures:
            N = len(structure)
            epsilon, distances, mask=self.get_matrices(structure, cut_radius, h)
            # epsilon = epsilon*mask
            r_hat = nu.zeros((N,N,3))
            for i in range(N):
                for j in range(N):
                    if not i == j:
                        vec = structure.positions[j] -structure.positions[i]
                        norm = nu.linalg.norm(vec)
                        r_hat[i,j] = vec/norm
            # if the given structure contains forces, use them, otherwise
            # assume that the structure is optimized
            try:
                forces_DFT = structure.get_forces()
            except:
                forces_DFT = nu.zeros((N,3))
            calc_new = self.solve_ground_state(structure, charge=charge, calc=calc)
            if not calc_new == None:
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
                        # FIXME this can be done easier/faster
                        indices = nu.array(mask[i,:]*epsilon[i,:])
                        indices -= 1
                        ppoints  = v_rep_points[indices]
                        f = forces_res[i].copy()
                        f -= nu.dot((mask[i,:]*ppoints),r_hat[i,:])
                        res += nu.linalg.norm(f)
                    return res

                v_rep_points, last_res_forces, minimized = self.find_forces(residual_forces, v_rep_points)
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
                    print "The minimization of forces did not converge!"
        if len(points) > 0:
            points = nu.array(points)
            # Normalize weights so that the best results corrensponds
            # to a weight given by the user (smallest residual force
            # after the minimization => largest weight)
            # +0.001 assures that 1/w does not diverge
            min_res_force = nu.min(points[:,2]) + 0.001
            if color == None:
                color = self.get_color()
            for data in points:
                self.append_point([data[0], data[1], weight*min_res_force/(data[2]+0.001), color, label], comment)
                label = '_nolegend_'


    def get_matrices(self, structure, cut_radius, h):
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


    def find_forces(self, function, v_rep_points):
        """
        Try to minimize the residual forces by finding matrix
        elements epsilon_ij = V'_rep(|r_ij|).
        """
        from scipy.optimize import fmin
        last_res_forces = 0
        N = len(v_rep_points)
        while True:
            print "Found %i different bond lengths." % N
            ret = fmin(function, v_rep_points, full_output=1)
            v_rep_points = ret[0]
            forces = ret[1]
            print "The sum of the norm of the net forces: %0.4f" % (forces)
            if ret[4] == 0 and nu.abs(forces-last_res_forces) < 0.001:
                return v_rep_points, forces, True
            last_res_forces = ret[1]
        return v_rep_points, last_res_forces, False


#    def append_energy_curve(traj,edft,charge,bonds,weight):
#        """ 
#        Fit V_rep'(r) into energy curve E(r) in given trajectory.
#        
#        Bonds are the atom index pairs (starting from zero) defining
#        the bonds that change (and should all be equal). Whatever the
#        system is, the ONLY missing energy contribution should thus
#        be the repulsive energy between given bonds! Beware of 
#        the dependency on other repulsive potentials (fitting should
#        depend only on the electronic part).
#        
#        UUSIKSI.
#        """
#        raise NotImplementedError('should (might) work; check!!')
#    
#        trajectory=Atoms(traj)
#        name=trajectory[0].get_name()
#        n=len(trajectory)
#        ebs=[]
#        nb=len(bonds)
#        rl=[]
#        # construct edft(r) function
#        for i,atoms in enumerate(trajectory):
#            ebs.append( self.get_energy(atoms,charge) )
#            lengths=[atoms.distance(b[0],b[1]) for b in bonds]
#            if max(lengths)-min(lengths)<1E-6:
#                print 'bond lengths:',lengths
#                raise ValueError('bond lengths vary too much.')
#            rl.append(nu.average(nu.array(lengths)))
#            
#        # now edft(r)=ebs(r)+nb*vrep(r); make vrep(r) and get vrep'(r)
#        vrep=[(edft[i]-ebs[i])/nb for i in range(n)]
#        v=SplineFunction(rl,vrep)
#        for r in rl:
#            self.add_point([r,v(r,der=1),weight/n,name],\
#                 comment='energy curve %s' %name)
#        raise NotImplementedError('should (might) work; check!!')


    def append_equilibrium_structure(atoms,relaxed=None):
        """
        Fit V_rep'(r) into given equilibrium (DFT) structure.
        
        atoms is a file name or Atoms object. The only requirement
        is that the ONLY missing force component missing
        from the TB calculation for the atoms in 'relaxed' list
        is the repulsion to be fitted. If relaxed is None, all
        the residual force is minimized for all atoms.
        
        TARKISTETTAVA.
        """
        from copy import copy
        raise NotImplementedError('still to be checked')
        if relaxed is None:
            relaxed=range(len(atoms))
        forces0=self.get_forces(atoms,forces=True)
        param=copy(self.param)
        #vrep=Spli...
        self.repulsion_forces(atoms,vrep)
        #residuals=get_residuals(atoms,forces0)
        #minimize residuals -> vrep'          
        #r
        #vrep(r,der=1)


    def write_fitting_data(self, filename):
        import pickle
        f = open(filename,'w')
        pickle.dump(self.deriv, f)
        f.close()
        f = open(filename[:-4]+'_human_readable'+filename[-4:],'w')
        print >> f, "# %0.6s %0.6s %0.6s %0.6s" % ('R', 'V', 'weight','info')
        for data in self.deriv:
            print >> f, "%0.6f %0.6f %0.3f # %s" % (data[0], data[1], data[2], data[4])
        f.close()


    def load_fitting_data(self, filename):
        import pickle
        f = open(filename,'r')
        self.deriv = pickle.load(f)
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
            structures.append(read(filename))
        else:
            raise Exception("Unknown file format: %s" % structure)
        return structures


    def get_color(self):
        color = self.colors[self.color_index]
        self.color_index += 1
        if self.color_index == len(self.colors):
            self.color_index = 0
        return color


class ParametrizationTest:
    """
    A tool to examine how well your parametrization agrees with
    given ase-trajectories.

    trajectories:  list of trajectories you want to compare
    charges:       the charges of the systems in trajectories
    """

    def __init__(self, pars, trajectories, charges=None):
        if charges != None:
            assert len(trajectories) == len(charges)
            self.charges = charges
        else:
            charges = nu.zeros(len(trajectories))
        self.trajectories = trajectories
        self.points = []
        self.ref_points = []
        self.pars = pars
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
        for i_traj, charge in zip(range(len(self.trajectories)), self.charges):
            pl.figure(i_traj)
            pl.title(self.trajectories[i_traj])
            if i_par == 0:
                self.plot_ref(i_traj)
            self.compare_trajectory(i_traj, charge, tables, i_par)


    def compare_trajectory(self, i_traj, charge, tables, i_par):
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
                atoms.set_calculator(Calculator(SCC=True, charge=charge, tables=tables))
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
        color = self.colors[i_par]
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


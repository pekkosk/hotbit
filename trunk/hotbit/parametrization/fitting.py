from hotbit import Element
from hotbit import Calculator
import numpy as nu
from ase.units import Bohr,Hartree
from box import mix
from box import Atoms
from box.interpolation import Function


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

        print 'r_dimer      =',r_dimer
        print '1.5 x r_dimer=',1.5*r_dimer


    def __call__(self,r,der=0):
        """ Return V_rep(r) or V_rep'(r) """
        if self.v is None:
            return self.parametric_potential(r,self.param,der=der)
        else:
            self.v(r,der=der)


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
        import pylab as pl
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
        pl.xlim(xmin=rmin)
        pl.ylim(ymin=self(rmin,der=1))
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
        
        
    def solve_ground_state(self, atoms, charge=None):
        """
        vahan kuten get_energy
                
        """
        from copy import copy
        calc = copy(self.calc)
        if charge != None:
            calc.set("charge",charge)
        atoms.set_calculator(calc)
        try:
            # FIXME make these output to file also
            atoms.get_potential_energy()
        except Exception:
            del(calc)
            calc = None
        return calc


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


    def fit2(self, r_cut):
        """
        Fit first derivative of 2nd order polynomial p(x)
        into points {r,V_rep'(r)}. The result is used as an initial
        guess in fitting the derivative of 3rd order p(x).
        Finally Nth order p(x) is given as an initial guess
        to the final fitting procedure.
        """
        self.r_cut = r_cut
        from scipy.optimize import fmin
        N = 5
        params = nu.zeros(2, float)
        for o in range(N-1):
            p = nu.zeros(len(params) + 1)
            p[:-1] = params
            params = p
            params = fmin(self.fitting_function2, params)
        self.param=nu.zeros(self.order,float)
        self.param[0:N+1] = params
        self.param = fmin(self.fitting_function,self.param,maxiter=self.maxiter, maxfun=self.maxiter)


    def fitting_function2(self,d):
        """
        Minimize this function in the fitting. The derivative of
        the polynomial should fit well to the points, it must be
        zero at r=r_cut and it should be relatively smooth curve.
        """
        ret = 0.0
        d[0] = 0
        d[1] = 0
        for point in self.deriv:
            r = point[0]
            dv = point[1]
            w = point[2]
            ret += ( dv-self.polynomial(r,d,der=1) )**2*w
        ret += self.curvature_penalty(d)
        return ret


    def polynomial(self, r, d, der=0):
        """
                                                 N                     
        Return the value of a polynomial p(r) = sum d[i]*(r_cut - r)**i
                                                i=0                    
        at point r.
        """
        r0 = self.r_cut
        if r>r0 or r<0:
            return 0.0
        if der == 0:
            return sum([d[i]*(r0-r)**i for i in range(len(d))])
        if der == 1:
            return sum([d[i]*(r0-r)**(i-1)*i*(-1) for i in range(1,len(d))])
        raise NotImplementedError("Only first derivative is available")


    def curvature_penalty(self, d):
        """
        Return a value between zero and +\infty that somehow describes
        the curvature of the polynomial (zero == straigth line)
        """
        rg = nu.linspace(0, self.r_cut, 100)
        ret = 0
        for i, r in enumerate(rg):
            if 1 <= i <= len(rg)-2:
                ret += ( (self.polynomial(rg[i+1],d) - 2*self.polynomial(r,d) + self.polynomial(rg[i-1],d)) / (rg[i+1]-rg[i])**2 )**2
        return ret/len(rg)


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
        weight:              Weight given for the points calculated
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
        import pylab
        from ase.io.trajectory import PickleTrajectory
        from ase import Atoms
        from copy import copy
        from box.interpolation import SplineFunction

        if not 'separating_distance' in kwargs:
            kwargs['separating_distance'] = 3.0
        if not 'h' in kwargs:
            kwargs['h'] = 1e-5
        if not 'weight' in kwargs:
            kwargs['weight'] = 1.0
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
        # sort the radii and corresponding energies to ascending order
        indices = R.argsort()
        R = R[indices]
        E_dft = E_dft[indices]
        E_bs = E_bs[indices]
        vrep = SplineFunction(R, (E_dft - E_bs)/N)

        if not 'color' in kwargs:
            color = 'b'
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
            self.append_point([r, vrep(r,der=1), kwargs['weight'], color, label], comment=comment)


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

from hotbit import Element
from hotbit import Calculator
import numpy as nu
from ase.units import Bohr,Hartree
from box import mix
from box import Atoms
from box.interpolation import Function


class RepulsiveFitting:
    
    def __init__(self,rep,r_cut=None,r_dimer=None,order=8,calc=None):
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
        self.param=nu.zeros(order,float)
        self.param[2]=10                    # initial guess...
        self.deriv=[]
        self.comments=''
        self.scale=1.025                    # scaling factor for scalable systems
        self.calc=calc
        self.v=None

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
            print>>o, r,self(r)
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
        pl.plot(r,vp,label='$V_{rep}(r)$')
        for s in self.deriv:
            #pl.scatter( [s[0]],[s[1]],s=100*s[2],label=s[3])
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
        
        
    def solve_ground_state(self, atoms, charge=0):
        """
        vahan kuten get_energy
                
        """
        from copy import copy
        calc = copy(self.calc)
        calc.set("charge",charge)
        atoms.set_calculator(calc)
        try:
            # FIXME make these output to file also
            atoms.get_potential_energy()
            return calc
        except Exception, ex:
            del(calc)
            raise Exception(ex)


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


    def energy_curve(self, dft_traj, elA, elB, **kwargs):
        """
        Calculates the V'rep(r) from a given DFT ase-trajectory for elements
        A and B:

                    E'_DFT(R) - E'_BS(R)
        V'rep(R) =  ------------------ ,
                            N

        where R is the distance between the elements in the system and
              N is the number of different A-B pairs that are taken
              into account.
        """
        # FIXME crashes if the dftb calculation does not converge
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
        traj = PickleTrajectory(dft_traj)
        R, E_dft, N = self.process_trajectory(traj, elA, elB, **kwargs)
        E_bs = nu.zeros(len(E_dft))
        M = 0
        if 'frames' in kwargs:
            frames = kwargs['frames']
        else:
            frames = len(E_dft)
        for i in range(frames):
            try:
                atoms=copy(traj[i])
                calc = self.solve_ground_state(atoms)
                E_bs[i] = calc.get_potential_energy(atoms)
                del(calc)
                M = i+1
            except Exception, ex:
                print ex
                print "Could not converge after %ith point." % M
                break
        traj.close()
        vrep = SplineFunction(R[:M], (E_dft[:M] - E_bs[:M])/N)

        if not 'color' in kwargs:
            color = 'b'
        else:
            color = kwargs['color']
        if not 'label' in kwargs:
            label = ''
        else:
            label = kwargs['label']
        for i, r in enumerate(R[:M]):
            if i > 0:
                label='_nolegend_'
            self.append_point([r, vrep(r,der=1), kwargs['weight'], color, label], comment="Point from energy curve fitting")
        if 'plot' in kwargs and kwargs['plot']:
            pylab.plot(R[:M], E_dft[:M], label='DFT')
            pylab.plot(R[:M], E_bs[:M], label='DFTB-Vrep')
            pylab.plot(R[:M], vrep(R[:M]), label="Fitted repulsion")
            pylab.plot(R[:M], vrep(R[:M], der=1), label="Derivative of the repulsion")
            pylab.plot(R[1:M-1], vrep(R[1:M-1],der=1), 'o', label="Added points")
            xmin, xmax, ymin, ymax = pylab.axis()
            pylab.plot((self.r_dimer, self.r_dimer),(ymin,ymax))
            pylab.legend()
            pylab.show()


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
        self.assert_fixed_bond_lengths_except(traj, elA, elB, **kwargs)
        for i, image in enumerate(traj):
            atoms = copy(Atoms(image))
            r, N = self.get_distance_of_elements(elA, elB, atoms, **kwargs)
            E_dft[i] = image.get_total_energy()
            R[i] = r
        indices = R.argsort()
        R = R[indices]
        E_dft = E_dft[indices]
        return R, E_dft, N


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
        for i in range(len(atoms)):
            for j in range(i,len(atoms)):
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
            for i in range(len(atoms)):
                for j in range(i,len(atoms)):
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
        for i in range(len(positions)):
            for j in range(i, len(positions)):
                a = positions[i]
                b = positions[j]
                if ( a.symbol == elA and b.symbol == elB ) or \
                   ( a.symbol == elB and b.symbol == elA ):
                    R.append(nu.linalg.norm(a.position - b.position))
        R.sort()
        R = nu.array(R)
        R_min = R[0]
        N = nu.sum(nu.where(nu.abs(R - R_min) < h, 1, 0))
        for i in range(N):
            if R[i] - R[0] > h:
                raise AssertionError("Element pairs have too much difference in their relative distances.")
        for r in R[N:]:
            if r < separating_distance:
                raise AssertionError("Element pairs are too close to each other.")
        return nu.average(R[:N]), N


    def fitting_function(self,d):
        """ Minimize this function in the fitting. """
        chi2=0.0
        for point in self.deriv:
            r=point[0]
            vp=point[1]
            w=point[2]
            chi2+=(vp-self.parametric_potential(r,d,der=1))**2*w
        return chi2
        

    def fit(self):
        """ Fit V_rep(r) into points {r,V_rep'(r)}. """
        from scipy.optimize import fmin
        self.param = fmin(self.fitting_function,self.param)
    
    
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
        self.append_scalable_system(dimer,0.0,weight,comment='dimer at %.4f %s' %(self.r_dimer,chr(197)))
        
        
    def append_scalable_system(self,system,charge,weight,comment=None):
        """ Use scalable equilibrium (DFT) system in repulsion fitting. 
        
        TARKISTA. (KAYTA VOIMIA SUORAAN?)
        """
        if type(system)==type(''):
            atoms=read(file)
            name=file
        else:
            atoms=system
            name=atoms.get_name()
        x=self.scale
        r=atoms.mean_bond_length()        
        bonds=atoms.number_of_bonds()        
        e1=self.get_energy(atoms,charge)
        atoms.scale_positions(x)
        e2=self.get_energy(atoms,charge)   
        
        dEdr=(e2-e1)/(x*r-r)
        self.add_point([r,-dEdr/bonds,weight,name])
        if comment is None:
            comment='scalable %s' %name
        self.add_fitting_comment(comment)                    
    
    
    def append_energy_curve(traj,edft,charge,bonds,weight):
        """ 
        Fit V_rep'(r) into energy curve E(r) in given trajectory.
        
        Bonds are the atom index pairs (starting from zero) defining
        the bonds that change (and should all be equal). Whatever the
        system is, the ONLY missing energy contribution should thus
        be the repulsive energy between given bonds! Beware of 
        the dependency on other repulsive potentials (fitting should
        depend only on the electronic part).
        
        UUSIKSI.
        """
        raise NotImplementedError('should (might) work; check!!')
    
        trajectory=Atoms(traj)
        name=trajectory[0].get_name()
        n=len(trajectory)
        ebs=[]
        nb=len(bonds)
        rl=[]
        # construct edft(r) function
        for i,atoms in enumerate(trajectory):
            ebs.append( self.get_energy(atoms,charge) )
            lengths=[atoms.distance(b[0],b[1]) for b in bonds]
            if max(lengths)-min(lengths)<1E-6:
                print 'bond lengths:',lengths
                raise ValueError('bond lengths vary too much.')
            rl.append(nu.average(nu.array(lengths)))
            
        # now edft(r)=ebs(r)+nb*vrep(r); make vrep(r) and get vrep'(r)
        vrep=[(edft[i]-ebs[i])/nb for i in range(n)]
        v=SplineFunction(rl,vrep)
        for r in rl:
            self.add_point([r,v(r,der=1),weight/n,name],\
                 comment='energy curve %s' %name)
        raise NotImplementedError('should (might) work; check!!')
    
    
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
        
    
 
    
    
    
        

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
        raise NotImplementedError('Not yet consistently implemented; however, you may comment this exception to play around.')
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
        """ Append Vrep into par-file, including comments. """
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
            pl.scatter( [s[0]],[s[1]],s=100*s[2],label=s[3])
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
        
        
    def get_energy(self,atoms,charge,forces=False):
        """
        Calculate energy (or forces) for given structure with given charge
        and given parameters (elements & tables as well).
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
        

    def fit(self):
        """ Fit V_rep(r) into points {r,V_rep'(r)}. """
        from scipy.optimize import fmin
        self.param = fmin(self.fitting_function,self.param)
    
    
    def repulsion_forces(self,atoms,vrep):
        """ 
        Return the repulsive forces for atoms using given
        vrep(r) function for present element pairs.
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
        """ Add point to vrep'-fitting: data=[r,v',w,info] """
        self.deriv.append(data)
        if comment is not None:
            self.add_fitting_comment(comment)           
           
           
    def append_dimer(self,weight): 
        """ Use dimer bond length in fitting. """
        dimer=Atoms(symbols=[self.sym1,self.sym2],\
                    positions=[(0,0,0),(self.r_dimer,0,0)],\
                    pbc=False,cell=[100,100,100])
        self.use_scalable_system(dimer,0.0,weight,comment='dimer at %.4f %s' %(self.r_dimer,chr(197)))
        
        
    def append_scalable_system(self,system,charge,weight,comment=None):
        """ Use scalable equilibrium (DFT) system in repulsion fitting. """
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
        
    
 
    
    
    
        
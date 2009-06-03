from hotbit import *
from hotbit import Calculator
import numpy as nu
from box import mix
from box import Atoms
import os
find=mix.find_value


class RepulsivePotential:
    def __init__(self,rep,r_cut=None,r_dimer=None,order=8,tables=None,elements=None):
        if type(rep)==type(''):
            raise NotImplementedError('read .par file')
        else:
            self.elm1=Element(rep[0])
            self.elm2=Element(rep[1])
            self.sym1=self.elm1.get_symbol()
            self.sym2=self.elm2.get_symbol()
            self.r_dimer=r_dimer    
            self.r_cut=r_cut                    # the cutoff radius
            self.r_small=r_dimer*0.5            # use ZBL repulsion for r<r_small
            self.order=order                    # order of Vrep polynomial 
            self.param=nu.zeros(order,float)
            self.param[2]=0.5                   # initial guess...
            self.deriv=[]
            self.comments=''
            self.scale=1.025                    # scaling factor for scalable systems
            
            self.elements={self.sym1:self.elm1.get_file(),self.sym2:self.elm2.get_file()}
            if elements!=None:
                self.elements.update(elements)
            self.tables={'%s%s' %(self.sym1,self.sym2):'%s_%s.par' %(self.sym1,self.sym2)}
            if tables!=None:
                self.tables.update(tables)
            self.hb={'SCC':True,'convergence':1E-9,'tables':self.tables,'elements':self.elements,\
                     'verbose':False,'txt':'fitting.cal'}
                     
            print 'r_dimer      =',r_dimer
            print '1.5 x r_dimer=',1.5*r_dimer
    def add_fitting_comment(self,s):
        self.fitting+='|'+s        
            
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
            
    def parametric_potential(r,d,der=0):
        """ 
        Return Vrep(r) with given parameter representation. 
        
        Vrep(r) =sum_i[0,order] d_i (r_c-r)**i
        Vrep'(r)=sum_i[1,order] d_i (r_c-r)**(i-1)*i*(-1)
        """
        re=0.0
        if r>self.r_cut or r<0:
            return 0.0
        else:
            if der==0:
                return sum( [d[i]*(self.r_cut-r)**i for i in range(self.order)] )
            elif der==1:
                return sum( [d[i]*(self.r_cut-r)**(i-1)*i*(-1) for i in range(1,self.order)] )                
               
    def write_to_par(self,txt=None,points=100,append=True):
        """ Write Vrep into .par file, including comments. """
        if txt is None:
            txt='%s_%s.par' %(self.elm1,self.elm2)
        mode=(append,'a','w')
        o=open(txt,mode)
        print>>o, 'fitting='
        print>>o, self.comments
        print>>o, '\n\nrepulsion='
        for r in nu.linspace(0.1,self.r_cut,points):
            print>>o, r,self(r)
        o.close()
    
    def use_dimer(self,weight):  
        dimer=Atoms(symbols=[self.sym1,self.sym2],positions=[(0,0,0),(self.r_dimer,0,0)],\
                    pbc=False,cell=[100,100,100])
        self.use_scalable_system(dimer,0.0,weight,comment='dimer at %7.4f %s' %(self.r_dimer,chr(197)))
        
    def use_scalable_system(self,system,charge,weight,comment=None):
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
        self.deriv.append([r,-dEdr/bonds,weight,name])
        if comment is None:
            comment='scalable %s' %name
        self.add_fitting_comment(comment)
            
    def scale_positions(self,x,cell=False):   
        """ Scale the whole system by x; also the unit cell. """
        self.set_cell(self.get_cell()*x,fix=False)
        self.reduce_atoms_into_cell()
        
    def get_energy(self,atoms,charge):
        """
        Calculate energy for given structure with given charge
        and given parameters (elements & tables as well).
        """
        if type(atoms)==type(''):
            atms=read(file)
        else:
            atms=atoms
        param=self.hb.copy()
        param.update({'charge':charge})
        calc=Calculator(**param)
        atms.set_calculator(calc)
        e=atms.get_potential_energy()
        calc.finalize()
        return e
       
    #def make_fitting(self):
        #rmin  = 0.5
        #rmax  = r_cut   
        #self.deriv.sort()
        #der = vec(self.deriv)
        
        #from scipy.optimize import fmin
        
        #func  = PFunction( rdV[:,0],rdV[:,1],rdV[:,2],r_cut ) 
        #param = fmin(func.chi2,param)
        
        #if Npar==1: param=[param]
        #print "Optimized parameters:",param
        #print "\nV_rep(r_cut-d)=sum_(i=1,N) -p_i*d**i; p_i=\n",param[0]**2,param[1:],'\n'
    
    
    
if __name__=='__main__':
    #elm=Element('H')
    #elm.write_to_file('koe')
    #print elm.__dict__
    rep=RepulsivePotential(['C','Au'],r_cut=2.6459,r_dimer=1.85899,tables={'AuC':'Au_C.par'})
    rep.use_dimer(1.0)
    #rep.write_to_file()
    print rep.__dict__
        
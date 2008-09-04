#!/usr/bin/env python
 
import simu
import numpy as nu
from box import mix
vector=nu.array

def simple_propagation(molecule,dt):
    """ Propagate molecule (given mass, velocity and forces) for dt. """
    for atom in molecule:
        r=atom.get_position()
        f=atom.get_force()        
        v=atom.get_velocity()
        m=atom.get_mass()
        r=r+v*dt+0.5*f/m*dt**2
        atom.set_position(r)


class FireMinimizer:
    """ Class for molecule optimization objects. """
    
    def __init__(self,calc,molecule,fcrit=5E-4):
        """ Optimize the given molecule. """
        self.molecule=molecule
        self.calc=calc
        self.N=self.molecule.get_N()
        self.dt=0.5
        self.it=0
        self.cuts=0
        self.last_cut=0
        self.fcrit=fcrit
        
        # FIRE parameters (standard)
        self.a_start=0.1
        self.a=self.a_start
        self.f_inc=1.1
        self.f_dec=0.5
        self.f_a=0.99
        self.N_min=4
        self.dt_max=1.0
        self.N_MIC=0
        self.v_mix=True
    
    def step(self):
        """ Take one MD minimization step. """
        self.it+=1
        verlet(1,self.molecule,self.dt,mass='unit')
        self.calc.update(self.molecule)
        verlet(2,self.molecule,self.dt,mass='unit')
        done=self.fire_step()
        return done
        
    
    def __call__(self):
        """ Make one MD time step. """
        return self.step()
        
    def fire_step(self):
        """ FIRE optimization propagation step. """
        if self.it<=self.N_MIC:
            for atom in self.molecule:
                fv=nu.vdot(atom.get_velocity(),atom.get_force())
                if fv<0: atom.set_velocity( vector([0.,0.,0.]) )
        else:
            f=self.molecule.get_forces(mode='global')
            v=self.molecule.get_velocities(mode='global')
            fv=nu.vdot(f,v)
            if self.v_mix: 
                v=(1-self.a)*v + self.a*f*mix.norm(v)/mix.norm(f)
            
            if fv<0.0:
                self.cuts+=1
                self.last_cut=self.it
                self.dt=self.dt*self.f_dec
                self.a=self.a_start
                self.molecule.set_velocities( vector([0.0]*self.N*3) )
            elif self.it-self.last_cut>self.N_min:
                self.dt=min( self.dt_max,self.dt*self.f_inc )
                self.a=self.a*self.f_a
        
        fmax=mix.max_norm( self.molecule.get_forces() )   
        if fmax<self.fcrit:
            return True
        else:
            return False             
            
    #def verlet(self,step):
            #""" The velocity-Verlet propagation step (1 or 2). """
            #for atom in self.molecule:
                #r=atom.get_position()
                #a=atom.get_force()    #MASS=1
                #v=atom.get_velocity()
                #dt=self.dt
                #if step==1:
                    #atom.set_position( r+v*dt+0.5*a*dt**2 )
                    #atom.set_velocity( v+0.5*a*dt )
                #elif step==2:
                    #atom.set_velocity( v+0.5*a*dt )

def verlet(step,molecule,dt,mass='default'):
        """ The velocity-Verlet propagation step (1 or 2). """
        for atom in molecule:
            r=atom.get_position()
            f=atom.get_force()    
            v=atom.get_velocity()
            if mass=='default':
                m=atom.get_mass()
            elif mass=='unit':
                m=1.0
            a=f/m
            if step==1:
                atom.set_position( r+v*dt+0.5*a*dt**2 )
                atom.set_velocity( v+0.5*a*dt )
            elif step==2:
                atom.set_velocity( v+0.5*a*dt )
                    
def pair_distribution_list(molecule):
    """ Return the array r_ij=|r_i-r_j| for all pairs a 1-D array. """
    r=molecule.get_positions()
    N=molecule.get_N()
    rij=[]
    for i in xrange(N):
        for j in xrange(i+1,N):
            rij.append( molecule.distance(molecule.get_atom[i],molecule.get_atom[j]) )
    return vector(rij)        
            
def pair_distribution_function(molecule,rmin,rmax,sigma=0.7):
    """
    Return the pair distribution function within [rmin,rmax] with Gaussian broadening.
    """
    rij=pair_distribution_list(molecule)
    rij=nu.sort(rij)
    grd=mix.grid(rmin,rmax,200)
    g  =vector([grd,nu.zeros(len(grd))]).transpose()
    for x in rij:
        for i in xrange(len(grd)):
            g[i,1]+=mix.gauss_fct(grd[i],mean=x,sigma=sigma)
    return g
        
def average_bond_length(molecule):
    """ Average bond length using pair distribution function. """
    g   = pair_distribution_function(molecule,0.0,15.0)
    max = False
    eps = 1e-6
    # find the first max and following min of pair distr. fct
    if molecule.get_N==2:
        r=molecule.get_positions()
        return molecule.distance(molecule.get_atom[0],molecule.get_atom[1])
    else:    
        for i in xrange(len(g)-1):
            if g[i,1]>g[i+1,1]: max=True
            if max==True and (g[i,1]<g[i+1,1] or g[i,1]<eps or i==len(g)-3):
                rcut = g[i,0]
                rij  = pair_distribution_list(molecule)
                return sum( nu.where(rij<=rcut,rij,0.0) )/sum( rij<=rcut )
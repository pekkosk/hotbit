import numpy as npy


class GeneralizedLennardJones:
    def __init__(self, epsilon=1.0, sigma=1.0):
        self.epsilon = epsilon
        self.sigma = sigma
        self.positions = None

    def update(self, atoms):
        #FIXME: pbc
        #assert not atoms.get_pbc().any()
        if (self.positions is None or
            (self.positions != atoms.get_positions()).any()):
            self.calculate(atoms)

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self._forces

    def get_stress(self, atoms):
        return npy.zeros((3, 3))
    
    def calculate(self, atoms):
        #TODO: more effective implementation! 
        """
        Force: Force per unit cell
        Energy: Energy per unit cell
        
        P.S. This is not the most effective implementation
        """
        #TODO: implement for not WedgeAtomsm
        #FIXME: test for not WedgeAtomsm
        #FIXME: test for not WedgeAtomsm
        
        n_range = atoms.get_number_of_cells()
        print 'FIXME: 3d! GenLJ_get_number_of_cells = ', n_range
        
        positions = atoms.get_positions()
        self.energy = 0.0
        self._forces = npy.zeros((len(atoms), 3))
        
        print "self._forces = ", self._forces
        
        precision = 1e-10
         
        #energy
        for i1, p1 in enumerate(positions):
            #FIXME: for forces
            #for i2, p2 in enumerate(positions[:(i1)]):
            #for i2, p2 in enumerate(positions[:(i1 + 1)]):
            for i2, p2 in enumerate(positions[:]):
                #FIXME: add up to 3 n's 
                F = 0.0
                for n0 in range(n_range[0]):
                    #--# if (atoms.symmetry_operation(i2, (n0, 0, 0)) != 
                    #--#    atoms.symmetry_operation(p2, (n0, 0, 0))).any():
                    #--#    raise Exception ("Error!")
                    
                    #p2_n = atoms.symmetry_operation(i2,(n0,0,0))
                    p2_n = atoms.symmetry_operation(p2, (n0, 0, 0))
                    #print 'p2_n = ', p2_n
                    #print 'p1 = ', p1
                    # V(R = 0) = 0
                    if (abs(p1 - p2_n) < precision).all():
                        continue
                    diff = p2_n - p1
                    
                    # Do NOT count diagonal terms twice  
                    #if i1 == i2:
                    #    e_factor = 0.5
                    #else: e_factor = 1
                    e_factor = 0.5
                    
                    d2 = npy.dot(diff, diff)
                    c6 = (self.sigma ** 2 / d2) ** 3
                    c12 = c6 ** 2
                    
                    #print "self.energy = ", self.energy
                    self.energy += 4 * self.epsilon * (c12 - c6) * e_factor
        #self.energy = self.energy / 2.0    
        
        # forces
        for i1, p1 in enumerate(positions):
            #F = 0.0
            for i2, p2 in enumerate(positions):
                for n0 in range(n_range[0]):
                    print "n0 =", n0
                    p2_n = atoms.symmetry_operation(p2, (n0, 0, 0))
                    #print 'p2_n = ', p2_n
                    #print 'p1 = ', p1
                    # V(R = 0) = 0
                    if (abs(p1 - p2_n) < precision).all():
                        continue
                    diff = p2_n - p1
                    
                    d2 = npy.dot(diff, diff)
                    c6 = (self.sigma ** 2 / d2) ** 3
                    c12 = c6 ** 2
                    
                    M1 = atoms.transform_der((n0, 0, 0));                       print "M1 = ", M1
                    #Matrix_n = (M1 + npy.identity(3)) / 2;                        print "Matrix_n = ", Matrix_n
                    ##Matrix_n = (M1 - npy.identity(3))/2;                        print "Matrix_n = ", Matrix_n
                    Matrix_n = (npy.identity(3));                        print "Matrix_n = ", Matrix_n
                    vector_n1 = npy.dot(Matrix_n.T, diff);                       print "vector_n1 = ", vector_n1
                    #vector_n2 = npy.dot(diff, Matrix_n);                       print "vector_n2 = ", vector_n2
                    if (npy.dot(Matrix_n.T, diff) - 
                        npy.dot(diff, Matrix_n) > precision).any():
                        raise Exceptiopn("Bug?")
                    dF = 24 * self.epsilon * (2 * c12 - c6) / d2 * vector_n1;    print "dF = ", dF
                    #F += dF; print "F = ", F
                    
                    # minus here (as for original J-L calc.); why?
                    self._forces[i1] = self._forces[i1] - dF ; 
                    #print self._forces[i1] 
                    # FIXME: is sign ok?
                    
            #self._forces[i1] = self._forces[i1] + F
            
            #self._forces[i1] = self._forces[i1] - F
            # FIXME: this ok?
            #if (i1!=i2):
            #self._forces[i2] += F
            print "self._forces = ", self._forces                  
            
            

        self.positions = positions.copy()
import numpy as np

from ase.optimize.optimize import Dynamics
from numpy.random import random

class contour_exploration(Dynamics):
    def __init__(self, atoms, 
                step_size = 0.5,
                parallel_drift = 0.2,
                energy_target = None,
                force_parallel_step_scale = None,
                remove_translation = False,
                use_FS = True,
                initialize_old = True, initialization_step_scale = 1e-2,
                use_target_shift = True, target_shift_previous_steps = 10,
                angle_limiting_max_step = None,
                seed = 60622,
                verbose = False, 
                trajectory=None, logfile=None,
                append_trajectory=False, loginterval=1):

        '''Perpendicular drift breaks orbits like dimer form, so they spin on new axes'''


        if force_parallel_step_scale is None:
            # a hureistic guess since most systems will overshoot when there is drift
            FPSS = 0.9 + 0.5*parallel_drift
        else:
            FPSS = force_parallel_step_scale


        self.verbose = verbose
        self.seed = seed
        self.remove_translation = remove_translation
        self.use_FS = use_FS
        #### for FS taylor expansion
        self.T = None
        self.Told = None
        
        self.N = None
        self.Nold = None
        
        self.r_old = None
        self.r     = None

        #########

        if energy_target == None:
            self.energy_target = atoms.get_potential_energy()
        else:
            self.energy_target = energy_target

        self.previous_energies = np.zeros(target_shift_previous_steps)

        # we need velocities or this won't run and will produce NaNs,
        # if none are provided we make random ones
        p = atoms.get_momenta()
        masses = atoms.get_masses()[:, np.newaxis]
        v = p/masses
        from ase import units
        if self.dot(v,v) < 1e-6:
            from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
            MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
            print("No Velocities found, random ones applied")


        if initialize_old:
            self.step_size     = step_size*initialization_step_scale
            self.max_step_size = step_size*initialization_step_scale
            self.parallel_drift = 0.0
            ## should this be 0.0 for a better initial curvature?
            ## Or at least smaller than 1.0? doesn't seem to matter much
            self.force_parallel_step_scale = 1.0 
            self.use_target_shift = False
            self.atoms = atoms
            self.step()

        
        self.step_size     = step_size
        self.max_step_size = step_size
        self.angle_limiting_max_step = angle_limiting_max_step
        self.parallel_drift = parallel_drift
        self.force_parallel_step_scale = FPSS
        self.use_target_shift = use_target_shift
        
        ## initizing the previous steps at the target energy slows 
        ## target shifting untill the system has had 
        ## 'baseline_shift_previous_steps' steps to equilibrate 
        ## this should prevent occilations
        self.previous_energies.fill(energy_target)
        
        #print(self.previous_energies)
        
        ## loginterval exists for the MolecularDynamics classe but not the more
        ## general Dynamics class
        Dynamics.__init__(self, atoms,
                                logfile, trajectory, #loginterval,
                                append_trajectory=append_trajectory,
                                )

    ## Required stuff for Dynamics
    def todict(self):
        return {'type': 'contour-exploration',
                'dyn-type': self.__class__.__name__,
                'stepsize': self.step_size}

    def irun(self, steps=50):
        """ Call Dynamics.irun and adjust max_steps """
        self.max_steps = steps + self.nsteps
        return Dynamics.irun(self)


    def run(self, steps=50):
        """ Call Dynamics.run and adjust max_steps """
        self.max_steps = steps + self.nsteps
        return Dynamics.run(self)


    def converged(self):
        """ MD is 'converged' when number of maximum steps is reached. """
        return self.nsteps >= self.max_steps


    ### CED specific stuff
    def dot(self,a,b):
        '''a dot product for 3 column vectors of atoms.get_forces'''
        return (a*b).sum()

    def vector_rejection(self,a,b):
        '''returns new vector that removes vector a's projection vector b'''
        aout = a - self.dot(a,b)/self.dot(b,b) *b
        return aout

    def unit_vect(self,a):
        return a / np.linalg.norm(a)

    def step(self, f=None):
        debug = False
        atoms = self.atoms

        if f is None:
            f = atoms.get_forces()

        self.r = atoms.get_positions()
        current_energy = atoms.get_potential_energy()

        ## update our history of self.previous_energies to include
        ## our current energy
        ## np.roll shifts the values to keep nice sequential ordering
        self.previous_energies = np.roll(self.previous_energies,1) 
        self.previous_energies[0] = current_energy
        if self.verbose: print('Previous Energies' , self.previous_energies)

        if self.use_target_shift:
            target_shift = self.energy_target - np.mean(self.previous_energies)
            # does this help balance the small steps never reaching? nope. makes std worse
            #target_shift = target_shift/self.force_parallel_step_scale
            if self.verbose: print('Energy Target [Shift]:',self.energy_target,'[%f]'% target_shift)
        else: 
            target_shift = 0.0
        deltaU = current_energy - ( self.energy_target + target_shift)

        ## get the velocity vector and old kinetic energy for momentum rescaling
        p = atoms.get_momenta()
        KEold = atoms.get_kinetic_energy()
        masses = atoms.get_masses()[:, np.newaxis]
        v = p/masses

        ### get force an and correction distance
        f_norm = np.linalg.norm(f)
        self.N = self.unit_vect(f)
        # maybe reduce parallel step component to enhance stability?
        delta_s_perpendicular =  (deltaU/f_norm)* self.force_parallel_step_scale # can be positive or negative
        
        # remove velocity  projection on forces
        w_parallel = self.vector_rejection(v,self.N)
        self.T = self.unit_vect(w_parallel)
        ### 

        if self.Nold is None or self.use_FS == False: 
            # we cannot apply the FS taylor expansion without a previous point
            # so we do a dumb jump to start
            w_perpendicular =  self.N * delta_s_perpendicular

            # create drift unit with no projection on N or T
            w_drift = random((len(atoms),3))
            w_drift = self.vector_rejection(w_drift,self.N)
            w_drift = self.vector_rejection(w_drift,self.T)
            w_drift = w_drift - w_drift.sum(axis = 0)/len(atoms) # removes wandering so systems don't wander
            D = self.unit_vect(w_drift)


            # Without the use of curvature there is no way to estimate the limiting step size
            self.step_size = self.max_step_size


            if abs(delta_s_perpendicular) < self.step_size:
                contour_step_size = np.sqrt(self.step_size**2 - delta_s_perpendicular**2 )

                drift_vector    = contour_step_size * self.parallel_drift * D
                parallel_vector = contour_step_size * np.sqrt(1- self.parallel_drift**2)*  self.T

                dr = parallel_vector + drift_vector + w_perpendicular
            else:
                dr = self.step_size/abs(delta_s_perpendicular) * w_perpendicular


        else:

             ####################  Frenet–Serret formulas
            # this should keep the dr clear of the constraints
            delta_r = self.r - self.rold
            delta_s = np.linalg.norm(delta_r) 
            # approximation of delta_s we use this incase an adaptive step_size algo get used 

            delta_T = self.T - self.Told
            delta_N = self.N - self.Nold
            dTds = delta_T/delta_s
            dNds = delta_N/delta_s
            #kappa_T = np.linalg.norm(dTds)
            kappa = np.linalg.norm(dNds) # normals are better since they are fixed to the reality of forces
            Nfs = dTds/kappa

            if self.angle_limiting_max_step is not None:
            
                phi = np.pi/180*self.angle_limiting_max_step
                self.step_size = np.sqrt(2-2*np.cos(phi))/kappa
                self.step_size = min(self.step_size, self.max_step_size)


            if debug:
                print('Told\n' , self.Told)
                print('T\n' , self.T)
                
                print('Nold\n', self.Nold)
                print('N\n', self.N)
                print('Nfs\n', Nfs)
            
            if self.verbose:
                print('Curvature, kappa' , kappa, 'Radius, 1/kappa', 1/kappa)
                print('T.Nfs', self.dot( self.T, Nfs), 'N.Nfs', self.dot(self.N,Nfs))
            
            
            if abs(delta_s_perpendicular) < self.step_size:
                contour_step_size = np.sqrt(self.step_size**2 - delta_s_perpendicular**2 )
                delta_s_parallel = np.sqrt(1- self.parallel_drift**2) * contour_step_size
                delta_s_drift    = contour_step_size*self.parallel_drift
                
                if self.verbose:
                    print('step_size %f in (perpendicular), [parallel], {drift}: (%f) [%f] {%f}' % \
                        (self.step_size, delta_s_perpendicular,delta_s_parallel, delta_s_drift))
                
                N_guess = self.N + dNds* delta_s_parallel
                T_guess = self.T + dTds* delta_s_parallel
                
                w_perpendicular = delta_s_perpendicular*( N_guess)
                
                w_parallel = delta_s_parallel * self.T * (1 - (delta_s_parallel * kappa)**2/6.0) \
                            + self.N * kappa/2 * delta_s_parallel**2
                #########
                drift = random((len(atoms),3))
                drift = self.vector_rejection(drift,N_guess)
                drift = self.vector_rejection(drift,T_guess)
                # removes net translation, so systems don't wander
                drift = drift - drift.sum(axis = 0)/len(atoms) 
                D = self.unit_vect(drift)
                w_drift = D*delta_s_drift
                #########

                if debug:
                    print('Tguess\n', T_guess)
                    print('Nguess\n', N_guess)
                    print('w_perpendicular\n',w_perpendicular)
                    print('w_parallel\n',w_parallel)
                    print('w_drift\n',w_drift)


                dr = w_perpendicular + w_parallel + w_drift
                dr = self.step_size * self.unit_vect(dr) 
                # because we guess our orthonormalization directions, 
                # we should renormalize to ensure a correct step size 
            else:
                w_perpendicular =  self.N * delta_s_perpendicular
                dr = self.step_size/abs(delta_s_perpendicular) * w_perpendicular
            
        ## now that dr is done, we check if there is translation
        if self.remove_translation:
            net_motion = dr.sum(axis = 0)/len(atoms)
            #print(net_motion)
            dr = dr - net_motion
            dr_unit = dr/np.linalg.norm(dr)
            dr = dr_unit*self.step_size



        ##  save old positions before update
        self.Told = self.T
        self.Nold = self.N
        self.rold = self.r

        ### the update section ###

        # if we have constraints then this will do the first part of the
        # RATTLE algorithm:
        atoms.set_positions(self.r + dr)
        p = (atoms.get_positions() - self.r) * masses #/ self.dt

        # We need to store the momenta on the atoms before calculating
        # the forces, as in a parallel Asap calculation atoms may
        # migrate during force calculations, and the momenta need to
        # migrate along with the atoms.
        atoms.set_momenta(p, apply_constraint=False)

        ## Now we get the new forces! 
        f = atoms.get_forces(md=True)
        

        ## rescaling momentum to maintain constant kinetic energy!
        pnew = atoms.get_momenta()
        KEnew = atoms.get_kinetic_energy()
        Ms = np.sqrt(KEold/KEnew) #Ms = Momentum_scale
        atoms.set_momenta(Ms*pnew)

        ### Normally this would be the second part of RATTLE will be done here like this:
        ### atoms.set_momenta(atoms.get_momenta() + 0.5 * self.dt * f)
        if self.verbose: print('')
        return f

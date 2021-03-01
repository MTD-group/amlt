import numpy as np
from numpy.random import default_rng
from ase.optimize.optimize import Dynamics
import time

class contour_exploration(Dynamics):
    def __init__(self, atoms, 
                maxstep = 0.5,
                parallel_drift = 0.2,
                energy_target = None,
                angle_limit = None,
                force_parallel_step_scale = None,
                remove_translation = False,
                use_FS = True,
                initialize_old = True, initialization_step_scale = 1e-2,
                use_target_shift = True, target_shift_previous_steps = 10,
                seed = 19460926, 
                verbose = False, 
                trajectory=None, logfile=None,
                use_tangent_curvature = False,
                force_consistent=None,
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
        self.force_consistent = force_consistent
        self.kappa = 0.0 # initializing so logging can work
        self.use_tangent_curvature = use_tangent_curvature
        #### for FS taylor expansion
        self.T = None
        self.Told = None
        
        self.N = None
        self.Nold = None
        
        self.r_old = None
        self.r     = None
        ###### initialize rng
        self.rng = default_rng(seed)

        #######
        if energy_target is None:
            self.energy_target = atoms.get_potential_energy(force_consistent = self.force_consistent)
        else:
            self.energy_target = energy_target

        # these need to be initialized before the initialize_old step so it doesn't crash
        self.previous_energies = np.zeros(target_shift_previous_steps)
        ## initizing the previous steps at the target energy slows 
        ## target shifting untill the system has had 
        ## 'baseline_shift_previous_steps' steps to equilibrate 
        ## this should prevent occilations
        self.previous_energies.fill(self.energy_target)


        # we need velocities or this won't run and will produce NaNs,
        # if none are provided we make random ones
        p = atoms.get_momenta()
        masses = atoms.get_masses()[:, np.newaxis]
        v = p/masses
        from ase import units
        if np.linalg.norm(v) < 1e-6:
            ### maxwell boltzman is fine for most cases, but unit testing, 
            ### we want deterministic results
            #from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
            #MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
            v = self.rand_vect((len(atoms),3)) # we have to pass dimension since atoms are not yet stored
            atoms.set_momenta(v/masses)
            print("No Velocities found, random velocities applied")


        if initialize_old:
            #self.step_size  = maxstep*initialization_step_scale
            self.maxstep    = maxstep*initialization_step_scale
            self.parallel_drift = 0.0
            ## should force_parallel_step_scale be 0.0 for a better initial curvature?
            ## Or at least smaller than 1.0? doesn't seem to matter much
            self.force_parallel_step_scale = 1.0 
            self.use_target_shift = False
            self.atoms = atoms
            self.step()
            

        
        #self.step_size = maxstep # this interfers with logging
        self.maxstep   = maxstep
        self.angle_limit = angle_limit
        self.parallel_drift = parallel_drift
        self.force_parallel_step_scale = FPSS
        self.use_target_shift = use_target_shift
        
        
        
        #print(self.previous_energies)
        
        ## loginterval exists for the MolecularDynamics class but not the more
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

    # do I need this for CED?
    def converged(self):
        """ MD is 'converged' when number of maximum steps is reached. """
        return self.nsteps >= self.max_steps


    def log(self ):

        #T = time.localtime() # time logging seems silly
        if self.logfile is not None:
            name = self.__class__.__name__
            if self.nsteps == 0:
                #args = (" " * len(name), "Step", "Energy_Target", "Energy", "Radius", "Step_Size", "Energy_Deviation_per_atom")
                args = ("Step", "Energy_Target", "Energy", "Curvature", "Step_Size", "Energy_Deviation_per_atom")
                msg = "# %4s %15s %15s %12s %12s %15s\n" % args
                self.logfile.write(msg)
            #args = (name, self.nsteps, T[3], T[4], T[5], e, ast, fmax)
            #msg = "%s:  %3d %02d:%02d:%02d %15.6f%1s %12.4f\n" % args
            
            e = self.atoms.get_potential_energy(force_consistent =  self.force_consistent )
            dev_per_atom = (e-self.energy_target)/len(self.atoms)
            args = ( self.nsteps, self.energy_target, e, self.kappa, self.step_size, dev_per_atom)
            msg = "%6d %15.6f %15.6f %12.6f %12.6f %24.9f\n" % args
            self.logfile.write(msg)

            self.logfile.flush()


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
        
    def rand_vect(self,dims=None):
        if dims is None:
            vect = self.rng.random((len(self.atoms), 3))
        else:
            vect =  self.rng.random(dims)
        return vect
        

    def step(self, f=None):
        debug = False
        atoms = self.atoms

        if f is None:
            f = atoms.get_forces()

        self.r = atoms.get_positions()
        current_energy = atoms.get_potential_energy(force_consistent = self.force_consistent)

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
            w_drift = self.rand_vect()
            w_drift = self.vector_rejection(w_drift,self.N)
            w_drift = self.vector_rejection(w_drift,self.T)
            w_drift = w_drift - w_drift.sum(axis = 0)/len(atoms) # removes wandering so systems don't wander
            D = self.unit_vect(w_drift)


            # Without the use of curvature there is no way to estimate the limiting step size
            self.step_size = self.maxstep


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
            if self.use_tangent_curvature:
                kappa = np.linalg.norm(dTds)
            else:
                # normals are better since they are fixed to the reality of forces
                # I see smaller forces and energy errors in bulk systems with normals
                kappa = np.linalg.norm(dNds) 
            self.kappa = kappa
            
            # on a perfect trajectory, the normal can be computed this way, 
            # I think the normal should always be tied to forces
            Nfs = dTds/kappa

            if self.angle_limit is not None:
            
                phi = np.pi/180*self.angle_limit
                self.step_size = np.sqrt(2-2*np.cos(phi))/kappa
                self.step_size = min(self.step_size, self.maxstep)

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
                drift = self.rand_vect()
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

        # rescale momentum to KEold to get the net momentum??

        # We need to store the momenta on the atoms before calculating
        # the forces, as in a parallel Asap calculation atoms may
        # migrate during force calculations, and the momenta need to
        # migrate along with the atoms.
        atoms.set_momenta(p, apply_constraint=False)

        ## Now we get the new forces! 
        f = atoms.get_forces(md=True)
        
        
        ####### I don't really know if removing md=True from above will break compatibility.
        f_constrained = atoms.get_forces()
        # but this projection needs the forces to be consistent with the constraints. 
        
        
        
        ## set new velocities perpendicular so they get logged properly in the trajectory files
        vnew = self.vector_rejection( atoms.get_momenta()/masses, f_constrained)
        #vnew = self.vector_rejection(atoms.get_velocities(), f)
        atoms.set_momenta(vnew*masses)
        
        ## rescaling momentum to maintain constant kinetic energy.
        KEnew = atoms.get_kinetic_energy()
        Ms = np.sqrt(KEold/KEnew) #Ms = Momentum_scale
        atoms.set_momenta(Ms*atoms.get_momenta())

        ### Normally this would be the second part of RATTLE will be done here like this:
        ### atoms.set_momenta(atoms.get_momenta() + 0.5 * self.dt * f)
        if self.verbose: print('')
        return f

import numpy as np

from ase.md.md import MolecularDynamics

from ase import units
from numpy.random import random

class contour_exploration(MolecularDynamics):
    def __init__(self, atoms, timestep=None,
                step_size = 0.5,
                perpendicular_drift = 0.2,
                parallel_step_scale = 1.8, 
                remove_translation = True,
                #max_parallel_fraction = 0.9,
                energy_target = None,
                seed = 60622,
                trajectory=None, logfile=None,
                loginterval=1, append_trajectory=False):

        '''Perpendicular drift breaks orbits like dimer form, so they spin on new axes
parallel_step_scale over corrects the shift to target energy which greatly helps acheiveing the target since all the contours are curved. at 1.8 it reduces the deviation from the target energy by about a factor of 2.'''



            
        self.step_size = step_size
        self.parallel_step_scale = parallel_step_scale
        self.perpendicular_drift = perpendicular_drift
        self.seed = seed
        self.remove_translation = remove_translation
        
        # we need velocities or this won't run and will produce NaNs,
        # if none are provided we make random ones
        p = atoms.get_momenta()
        masses = atoms.get_masses()[:, np.newaxis]
        v = p/masses
        if self.dot(v,v) < 1e-6:
            from ase import units
            from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
            MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
            print("No Velocities found, random ones applied")

        from ase import units
        timestep = 1.0 * units.fs # not used but MD needs it
        MolecularDynamics.__init__(self, atoms, timestep, trajectory, logfile,
                                   loginterval,
                                   append_trajectory=append_trajectory,
                                   )


        if energy_target == None:
            self.energy_target = atoms.get_potential_energy()
        else:
            self.energy_target = energy_target


    def dot(self,a,b):
        '''a dot product for 3 column vectors of atoms.get_forces'''
        return (a*b).sum()
        
    def vector_rejection(self,a,b):
        '''returns new vector that removes vector a's projection vector b'''
        aout = a - self.dot(a,b)/self.dot(b,b) *b
        return aout


    def step(self, f=None):

        atoms = self.atoms

        if f is None:
            f = atoms.get_forces()



            

        deltaU = atoms.get_potential_energy() - self.energy_target

        p = atoms.get_momenta()
        masses = atoms.get_masses()[:, np.newaxis]
        v = p/masses
        
        f_over_norm_sqr = f / (np.linalg.norm(f)**2)
        
        # maybe reduce parallel step component to enhance stability?
        w_parallel = deltaU*f_over_norm_sqr * self.parallel_step_scale 
        
        #print('v',v)
        #print('f',f)
        
        
        #w_perpendicular = v - (v*f).sum()*f_over_norm_sqr
        w_perpendicular = self.vector_rejection(v,f)
        ###
        
        w_drift = random((len(atoms),3))
        w_drift = self.vector_rejection(w_drift,f)
        w_drift = self.vector_rejection(w_drift,w_perpendicular)
        w_drift = w_drift - w_drift.sum(axis = 0)/len(atoms) # removes wandering so systems don't wander
        
        w_drift_unit = w_drift / np.linalg.norm(w_drift) 
        
        
        #print(self.dot(w_drift,w_perpendicular), self.dot(w_drift,f) )
        
        w_parallel_norm = np.linalg.norm(w_parallel)
        w_perpendicular_unit = w_perpendicular/np.linalg.norm(w_perpendicular)
        
        if w_parallel_norm < self.step_size:
            perpendicular_step_size = np.sqrt(self.step_size**2 - w_parallel_norm**2 )
            
            drift_vector         = perpendicular_step_size * self.perpendicular_drift * w_drift_unit
            perpendicular_vector = perpendicular_step_size * np.sqrt(1- self.perpendicular_drift**2)*  w_perpendicular_unit
            
            dr = perpendicular_vector + drift_vector+ w_parallel
        else:
            dr = self.step_size/w_parallel_norm * w_parallel
        
        
        if self.remove_translation:
            net_motion = dr.sum(axis = 0)/len(atoms)
            #print(net_motion)
            dr = dr - net_motion
            dr_unit = dr/np.linalg.norm(dr)
            dr = dr_unit*self.step_size
        
        
        r = atoms.get_positions()

        # if we have constraints then this will do the first part of the
        # RATTLE algorithm:
        atoms.set_positions(r + dr)
        #if atoms.constraints:
        #    p = (atoms.get_positions() - r) * masses / self.dt
            
        
        p = (atoms.get_positions() - r) * masses #/ self.dt

        # We need to store the momenta on the atoms before calculating
        # the forces, as in a parallel Asap calculation atoms may
        # migrate during force calculations, and the momenta need to
        # migrate along with the atoms.
        atoms.set_momenta(p, apply_constraint=False)

        f = atoms.get_forces(md=True)
        

        # Second part of RATTLE will be done here:
        #atoms.set_momenta(atoms.get_momenta() + 0.5 * self.dt * f)
        return f





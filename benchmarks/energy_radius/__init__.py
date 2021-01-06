
from ase import io
import numpy as np

import os

def read_energy_radius_traj(file_name):

        traj = io.Trajectory(os.path.abspath(file_name),'r')
        #data = [ (im.get_volume()/len(im), im.get_potential_energy(force_consistent = True)/len(im)) for im in traj]
        data = [ (np.linalg.norm(im[0].position - im[1].position),
                im.get_potential_energy()/len(im)) for im in traj]
        traj.close()
        data = np.asarray(data).T
        smap = np.argsort(data[0])

        return np.array([data[0][smap],data[1][smap]] )
        
        
def read_force_radius_traj(file_name):

        traj = io.Trajectory(os.path.abspath(file_name),'r')
        #data = [ (im.get_volume()/len(im), im.get_potential_energy(force_consistent = True)/len(im)) for im in traj]
        data = [ (np.linalg.norm(im[0].position - im[1].position),
                np.linalg.norm(im.get_forces()[0]) ) for im in traj]
        traj.close()
        data = np.asarray(data).T
        smap = np.argsort(data[0])

        return np.array([data[0][smap],data[1][smap]] )

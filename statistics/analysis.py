
import numpy as np
from ase import io


def cancel_net_force(atoms):
    forces = atoms.calc.results['forces']
    drift = np.sum(forces,  axis = 0)/len(atoms)
    atoms.calc.results['forces'] = forces-drift



def read_evaluation_data(filename,
    MLIP_cache_dir = 'MLIP_cache',
    struct_types = [ 'known', 'polymorphD3', 'random' ],
    dyn_types    = [ 'md','relax','sp','ce','dimer'],
    use_forces = True,
    remove_force_drift = True,
    bad_data_traj_list = [],
    ):
    
    
    with open(filename,'r') as fid:
        lines = fid.readlines()
    
    image_pairs = []
    
    for line in lines:
        if '#' not in line:
            
            sline = line.split()
            struct_type = sline[-3]
            dyn_type    = sline[-2]
            
            if dyn_type in dyn_types and struct_type in struct_types:
                data_file = sline[-1]
                im_index = int(sline[-4])
                if data_file in bad_data_traj_list:
                    print('%s[%i] in bad_data_traj_list, skipping...' %(data_file, im_index))
                else:
                    
                    dirlist = data_file.split('/') 
                    cache_path = MLIP_cache_dir + '/' + dirlist[-4] + '/' + dirlist[-3] + '/'+ dirlist[-2] + '/' 
                    cache_file = cache_path+'%i.traj'%im_index
                    
                    data_traj = io.Trajectory(data_file,'r')
                    cache_traj = io.Trajectory(cache_file,'r')
                    
                    data_im = data_traj[im_index]
                    cache_im = cache_traj[0]
                    
                    if use_forces:
                        if remove_force_drift:
                            cancel_net_force(cache_im)
                            cancel_net_force(data_im)
                    
                    image_pairs.append( (cache_im, data_im, cache_path ) )

                    data_traj.close()
                    cache_traj.close()
    
    return image_pairs

#### #

def get_energy_lists(image_pairs):
    cache_energy = np.zeros(len(image_pairs))
    data_energy = np.zeros(len(image_pairs))
    for index, image_pair in enumerate(image_pairs):
        cache_energy[index] = image_pair[0].get_potential_energy()/len(image_pair[0])
        data_energy[index]  = image_pair[1].get_potential_energy()/len(image_pair[1])
    
    return cache_energy, data_energy

def get_force_list(image_pairs):
    cache_forces = []
    data_forces  = []
    for index, image_pair in enumerate(image_pairs):
        cache_forces.append( image_pair[0].get_forces())
        data_forces.append(  image_pair[1].get_forces())
    
    return cache_forces, data_forces

def get_composition_and_atom_count_arrays(image_pairs, elements):
    composition_array = np.zeros((len(image_pairs), len(elements)))
    atom_counts =       np.zeros(len(image_pairs))
    for index, image_pair in enumerate(image_pairs):
        atom_counts[index] = len(image_pair[1])
        symbols = image_pair[1].get_chemical_symbols()
        for el_index, el in enumerate(elements):
             composition_array[index,el_index] = symbols.count(el)/float(atom_counts[index])
    
    return composition_array, atom_counts
        
#############
def compute_force_norms_by_atom(forces):
    atom_force_norms = []
    for index in range(len(forces)):
        atom_force_norm = np.sqrt(   (forces[index]**2) .sum(axis=1) )
        atom_force_norms.append(atom_force_norm)
    return atom_force_norms


def compute_force_norms_by_image(forces):
    force_norms = []
    for index in range(len(forces)):
        force_norm = np.sqrt(  np.mean( (forces[index]**2) .sum(axis=1) ) )
        force_norms.append(force_norm)
    return np.array(force_norms)

## these all operate on the cache_forces and data_forces



def compute_relative_force_error_list(cache_forces, data_forces):
    relative_force_error_list = []
    for index in range(len(cache_forces)):
        relative_force_error_list.append((cache_forces[index] - data_forces[index])/(np.absolute(data_forces[index] )+0.01))
    return relative_force_error_list

def compute_rms_relative_force_error_by_atom(cache_forces, data_forces):
    relative_force_error_list = compute_relative_force_error_list(cache_forces, data_forces)
    rms_relative_force_error_by_atom = []
    for index in range(len(relative_force_error_list )):
        atom_force_error = np.sqrt(   (relative_force_error_list [index]**2) .sum(axis=1) )
        rms_relative_force_error_by_atom.append(atom_force_error)
    return rms_relative_force_error_by_atom


def compute_rms_relative_force_error(cache_forces, data_forces):
    rms_relative_force_error_by_atom = compute_rms_relative_force_error_by_atom(cache_forces, data_forces)
    rms_relative_force_error_list = []
    for index in range(len(rms_relative_force_error_by_atom )):
        rms_relative_force_error_list.append( np.sqrt(np.mean(rms_relative_force_error_by_atom[index]**2)) )
    return rms_relative_force_error_list
    
################
def compute_force_error_list(cache_forces, data_forces):
    force_error_list = []
    for index in range(len(cache_forces)):
        force_error_list.append(cache_forces[index] - data_forces[index])
    return force_error_list

    
def compute_rms_force_error_by_atom(cache_forces, data_forces):
    force_error_list = compute_force_error_list(cache_forces, data_forces)
    rms_force_error_by_atom = []
    for index in range(len(force_error_list)):
        atom_force_error = np.sqrt(   (force_error_list[index]**2) .sum(axis=1) )
        rms_force_error_by_atom.append(atom_force_error)
    return rms_force_error_by_atom

def compute_rms_force_error_by_image(cache_forces, data_forces):
    rms_force_error_by_atom = compute_rms_force_error_by_atom(cache_forces, data_forces)
    rms_force_error_list = []
    for index in range(len(rms_force_error_by_atom )):
        rms_force_error_list.append( np.sqrt(np.mean(rms_force_error_by_atom[index]**2)) )
    return np.array(rms_force_error_list)
################## force cosines

def compute_force_cosines_by_atom(cache_forces, data_forces):
    force_cosines_by_atom = []
    delta = 1e-6
    for index in range(len(cache_forces)):
        fmags = np.linalg.norm(cache_forces[index] ,axis = 1) * \
            np.linalg.norm(data_forces[index],axis = 1)
    
        fdot = (cache_forces[index]*data_forces[index]).sum(axis=1)
        
        force_cosines = fdot/(fmags+delta) 
        force_cosines_by_atom.append(force_cosines)
    return force_cosines_by_atom

def compute_force_cosines_by_image(cache_forces, data_forces):
    force_cosines_by_image = []
    for index in range(len(cache_forces)):
        fmags = numpy.linalg.norm(cache_forces[index]) * \
            numpy.linalg.norm(data_forces[index])
    
        force_cosines = (cache_forces[index]*data_forces[index]).sum()/fmags 
        force_cosines_by_image.append(force_cosines)
    return np.array(force_cosines_by_image)
#########################
def collapse_sub_lists(list_of_arrays):
    flat_list = []
    for ar in list_of_arrays:
        flat_list.append(ar.flatten())
    
    return np.concatenate(flat_list)

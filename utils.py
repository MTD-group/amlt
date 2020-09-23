

def reorder_image_list_for_balanced_atom_counts(image_list, ncores = 4):

    image_list.sort(key=lambda image: len(image), reverse=True)
    i = 0
    cores_and_atoms_list = []

    while i < len(image_list) and i < ncores: # the first pass fills each core
        core_index = i
        cores_and_atoms_list.append( [core_index, len(image_list[core_index]) ] )
        image_list[i].core_index = core_index
        i += 1 
    

    while i < len(image_list):
        cores_and_atoms_list.sort(key=lambda core_and_atoms: core_and_atoms[1] )
        core_index = cores_and_atoms_list[0][0] # the first one has the fewest atoms
        image_list[i].core_index = core_index
        ## update atoms 
        cores_and_atoms_list[0][1] += len(image_list[i])
        i += 1

    #print(cores_and_atoms_list)
    image_list.sort(key=lambda im: im.core_index)


    cores_and_atoms_list.sort(key=lambda core_and_atoms: core_and_atoms[0] )
    #print(cores_and_atoms_list)     

    atom_counts = [0 for i in range(ncores) ] # zeros
    for image in image_list:
        atom_counts[image.core_index] += len(image)

    return atom_counts
    
    




def get_image_list( basedirs = [''], 
    image_skip = 2, image_skip_offset = 0,
    traj_skip =1, traj_skip_offset = 0,
    image_offset = 0,
    traj_offset  = 0,  
    struct_types = [ 'random',   'known', 'polymorphD3'] ,
    dyn_types    = [ 'md', 'relax', 'sp', 'ce', 'dimer'] ,
    max_energy_per_atom = None, max_force_on_atom = None,
    remove_force_drift_in_training_data = True, return_file_paths = False):

    from ase import io
    from glob import glob
    import time
    import os
    import numpy as np

    
       
    def composition_str(sorted_elements, counts):
        comp_format = '( '
        for el, cnt in zip(sorted_elements, counts):
            comp_format = comp_format + el + '_%i '%cnt
        comp_format = comp_format + ')'

        return comp_format


    def remove_force_drift(atoms):
        forces = atoms.calc.results['forces']
        drift = np.sum(forces,  axis = 0)/len(atoms)
        atoms.calc.results['forces'] = forces-drift
        #return drift


    element_set = set()

    image_list = []
    file_path_list = []
    total = 0
    time1 = time.time()
    for basedir in basedirs:
        for struct_type in struct_types:
            for dyn_type in dyn_types:    
                top_direct = os.path.abspath( basedir)+ ('/%s_%s/')%(struct_type, dyn_type)
                #print(top_direct)
                if os.path.isdir(top_direct):
                        print(top_direct)
                        sub_total = 0
                        sub_dirs = sorted(glob(top_direct+'*/'))
                        #file_list.sort()
                        sub_dirs.sort(key= lambda x: len(x))
                        for sub_dir in sub_dirs:
                            name = sub_dir.split('/')[-2]
                            if name.isdigit():
                            
                                if int(name) >= traj_offset and int(name)%traj_skip == traj_skip_offset:
                                    traj = io.Trajectory(filename = sub_dir  + 'images.traj', mode='r')
                                    subsub_total = 0
                                    
                                    print (sub_dir.ljust(22), end = '')
                                    # for printing the compositions
                                    symbols = traj[0].get_chemical_symbols()
                                    element_set.update(symbols)
                                    sorted_elements = sorted(list(element_set))
                                    comp = [ symbols.count(el) for el in sorted_elements]
                                    comp_format = composition_str(sorted_elements, comp)
                                    
                                    for image_index in range(len(traj)):
                                        image = traj[image_index]
                                        
                                        #if image.get_potential_energy()/len(image) <= max_energy_per_atom:

                                        training_image = False

                                        if struct_type == 'known': training_image = True    
                                        #if (int(name)%trajskip) == trajskip_offset:
                                        if image_index >= image_offset and ((image_index) % image_skip) == image_skip_offset:
                                            training_image = True

                                        if training_image: 
                                            ## now we test for unreasonably high energies/forces 
                                            ## if thought we'd include the structure
                                            if max_energy_per_atom is not None:
                                                energy_per_atom = image.get_potential_energy()/len(image)
                                                if  energy_per_atom > max_energy_per_atom:
                                                    training_image = False
                                                    print('image', image_index, 'Energy too high:', energy_per_atom )
                                                    
                                            if max_force_on_atom is not None:
                                                max_force = np.linalg.norm(image.get_forces(), axis = 1).max()
                                                if max_force > max_force_on_atom:
                                                    training_image = False
                                                    print('image', image_index, 'Max Force too high:', max_force)
                                        #if (int(name) in bad_polymorphs) and struct_type == 'polymorphD3' :
                                        #    training_image = False

                                        if training_image:
                                            #print ('adding image no', image_index)
                                            if remove_force_drift_in_training_data: remove_force_drift(image)
                                            image_list.append(image)
                                            file_path_list.append([image_index, struct_type, dyn_type, sub_dir  + 'images.traj'])
                                            subsub_total += 1
                                            sub_total    += 1
                                            total        += 1
                                    
                                    print( (' atoms %i'%len(traj[0])).ljust(12) +\
                                            ('images loaded %i/%i '%(subsub_total, len(traj))).ljust(24) +comp_format )

                                
                                    traj.close()
                        # subtotal for this struct+dyn_type
                        print('sub_total: %i \n'% sub_total)
                        



    sum_total_atoms = 0
    for image in image_list:
        sum_total_atoms += len(image)


    print('Total Number of Images:', total)
    print('Total Atoms: %i' % sum_total_atoms)
    print('Time for file parsing is: {:.3f} sec.'.format(time.time() - time1))

    if return_file_paths:
        return image_list, file_path_list
    else:
        return image_list
        

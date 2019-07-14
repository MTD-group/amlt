"""Commands for generating DFT input files."""
import math
import numpy as np
import os
from amlt import reasonable_random_structure_maker, polymorphD3, try_mkdir
from ase import io
from os.path import isfile



def vasp_job_maker(name_prefix, 
                   jobs, 
                   job_command, 
                   job_script_name, 
                   job_script_template, 
                   md_temperature_range=(100, 600), 
                   submit = False,
                   random_structure_parameters = None, 
                   known_structures=[], 
                   polymorphD3_parameters = None,
                   first_structure = 'POSCAR.initial', 
                   magmom_filename = 'MAGMOMS.initial'):
    """Creates VASP input files and controls job submission.
    Args:
        name_prefix (String):
        jobs (Iterable):
        job_command:
        job_script_name:
        job_script_template:
        md_temperature_range (Tuple of numbers):
        submit (Boolean): Send jobs to SLURM scheduler if True
        random_structure_parameters: 
        known_structures (List): 
        polymorphD3_parameters:
        first_structure (String): Saved filename of starting structure
        magmom_filename (String): Saved filename of starting magnetic moments
    """

    
    twd = os.getcwd()

    for job_type in jobs:

        job_type_dir = job_type[1] +'_' + job_type[2] +'/'

        try_mkdir(job_type_dir)
        
        if job_type[1] =='known':
            n_structures = len(known_structures)
        else:
            n_structures = job_type[0]

        for structure_number in range(n_structures):

            struct_dir =job_type_dir + str(structure_number) +'/'
            try_mkdir(struct_dir)


            if not isfile(struct_dir+ first_structure):


                if job_type[1] == 'random':
                    atoms = reasonable_random_structure_maker(**random_structure_parameters)
                    #if callable(random_structure_parameters['magmom_generator']):
                    
                elif job_type[1] == 'polymorphD3':
                    index = np.random.random_integers(len(known_structures)-1)
                    atoms = polymorphD3(known_structures[index], **polymorphD3_parameters)
                    
                elif job_type[1] =='known':
                    atoms = known_structures[structure_number]
                else: 
                    raise Exception('Structure type "%s" not recognized'%job_type[1] )
                
                io.write(struct_dir+ first_structure, atoms, format = 'vasp')
                magmoms = atoms.get_initial_magnetic_moments()
                # magmom check
                if np.sqrt(magmoms.dot(magmoms)) > 0.0000001:
                    np.savetxt(struct_dir + magmom_filename,  magmoms.T)
                
                if job_type[2] == 'md':
                    temp = np.random.rand()*(max(md_temperature_range) - min(md_temperature_range)) + min(md_temperature_range)
                    np.savetxt(struct_dir+'temperature.txt', [temp,temp])
                    
                print(struct_dir, 'structure created')
                
                

            if not isfile(struct_dir+'OUTCAR'):

                fid = open(struct_dir+ job_script_name,'w')
                job_name = "{}_{}_{}_{}".format(
                        name_prefix,
                        job_type[1],
                        job_type[2],  
                        structure_number)
                fid.write(job_script_template.format(job_name, job_type[2])) 
                fid.close()

                    
                if submit:
                    os.chdir(struct_dir)
                    os.system('sbatch '+job_script_name)
                    os.chdir(twd)
                    print(struct_dir, 'job submitted')
                else:
                    print(struct_dir, 'job not yet run')
                
            else:
                images = io.read(struct_dir+'OUTCAR', index = ':')
                my_traj = io.trajectory.Trajectory( struct_dir + 'images.traj', mode = 'w')

                for atoms in images:
                    my_traj.write(atoms = atoms)
                my_traj.close()

                print(struct_dir.ljust(25), 'done with %i images'% len(images))







def kgrid_from_cell_volume(atoms, kpoint_density ):
    kpd = kpoint_density
    lengths_angles = atoms.get_cell_lengths_and_angles()
    #if math.fabs((math.floor(kppa ** (1 / 3) + 0.5)) ** 3 - kppa) < 1:
    #    kppa += kppa * 0.01
    lengths = lengths_angles[0:3]
    ngrid = kpd/atoms.get_volume() # BZ volume = 1/cell volume (withot 2pi factors)
    mult = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1 / 3)

    num_divf = [int(math.floor(max(mult / l, 1))) for l in lengths]
    kpdf = atoms.get_volume()*num_divf[0]*num_divf[1]*num_divf[2]
    errorf = abs(kpd - kpdf)/kpdf # not a type being 0.5 is much worse than being 1.5


    num_divc = [int(math.ceil(mult / l)) for l in lengths]
    kpdc = atoms.get_volume()*num_divc[0]*num_divc[1]*num_divc[2]
    errorc = abs(kpd - kpdc)/kpdc #same

    if errorc < errorf :
        num_div = num_divc
    else:
        num_div = num_divf

    #num_div = num_divf
    return num_div

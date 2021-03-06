"""Commands for generating DFT input files. Not intended to be called from
places besides control_script.py."""
import numpy as np
import os
from amlt import reasonable_random_structure_maker, PolymorphD3, try_mkdir
from ase import io
from os.path import isfile



def convert_to_traj(filename, traj_name= 'images.traj'):
    try:
        if '.traj' in filename:
            images = io.trajectory.Trajectory( filename, mode = 'r')
        else:
            images = io.read(filename, index = ':')
    
        output_traj = io.trajectory.Trajectory( traj_name, mode = 'w')
        for atoms in images:
            output_traj.write(atoms = atoms)
        output_traj.close()
    
    except:
        print('conversion failed! File may be incomplete.')
        images = []

    return images



def outcar_to_traj(outcar_name = 'OUTCAR', traj_name='outcar.traj'):
    return convert_to_traj(filename = outcar_name, traj_name= traj_name)




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
                   magmom_filename = 'MAGMOMS.initial',
                   rebuild_traj_cache = False):
    """Creates VASP input files and controls job submission.
    Args:
        name_prefix (String):
        jobs (Iterable of Lists): Lists contain
            [0]: the number of structures as an integer,
            [1]: the type of structural modification as a string,
                options are "known" which does not change the known crystalline
                polymorphs, "polymorphD3" which adds vacancies and atomic
                displacements to the known crystalline polymorphs, and "random"
                which creates totally random arrangements of atoms.
            [2]: which ase calculator template to use as a string,
                options are "sp", "md", or "relax". sp means single point
        job_command (String): command used to submit jobs to scheduler
        job_script_name (String): path to write the job script
        job_script_template (String): Text of job submission script with string
            formatting for calling the right calculator script.
        md_temperature_range (Tuple of numbers): Min and max MD temperature
        submit (Boolean): Send jobs to SLURM scheduler if True
        random_structure_parameters (Dict): Control parameters for making
            random structures.  See rrsm.py for details
        known_structures (List): List of CIF filepaths for known polymorphs
        polymorphD3_parameters (Dict): Control parameters for creating
            distorted structures with vacancies and displacements.
            See polymorphD3.py for more details
        first_structure (String): Saved filename of starting structure
        magmom_filename (String): Saved filename of starting magnetic moments
    """


    twd = os.getcwd()

    total_images = 0
    for job_type in jobs:

        job_type_dir = job_type[1] +'_' + job_type[2] +'/'

        try_mkdir(job_type_dir)

        if job_type[1] =='known':
            n_structures = len(known_structures)
        else:
            n_structures = job_type[0]
        
        # set outputfile that will be used for checking job progress
        if len(job_type)>=4:
            outputfile = job_type[3]
        else: #default to outcar if nto given
            outputfile = 'OUTCAR'

        # now we can loop over structures
        job_type_total_images = 0
        for structure_number in range(n_structures):

            struct_dir =job_type_dir + str(structure_number) +'/'
            try_mkdir(struct_dir)

            # If input structure files have not been generated, we need to
            # write a POSCAR and MAGMOM file based on the job type
            if not isfile(struct_dir+ first_structure):


                if job_type[1] == 'random':
                    atoms = reasonable_random_structure_maker(**random_structure_parameters)
                    #if callable(random_structure_parameters['magmom_generator']):

                elif job_type[1] == 'polymorphD3':
                    assert len(known_structures)>0
                    index = np.random.random_integers(low = 0, high=len(known_structures)-1)
                    atoms = PolymorphD3(known_structures[index], **polymorphD3_parameters).atoms_out

                elif job_type[1] =='known':
                    atoms = known_structures[structure_number]
                else:
                    raise Exception('Structure type "{}" not recognized'.format(job_type[1]))

                io.write(struct_dir+ first_structure, atoms, format = 'vasp', vasp5=True)
                magmoms = atoms.get_initial_magnetic_moments()
                # magmom check
                if np.sqrt(magmoms.dot(magmoms)) > 0.0000001:
                    np.savetxt(struct_dir + magmom_filename,  magmoms.T)

                #if job_type[2] == 'md':
                #    temp = np.random.rand()*(max(md_temperature_range) - min(md_temperature_range)) + min(md_temperature_range)
                #    np.savetxt(struct_dir+'temperature.txt', [temp,temp])

                print(struct_dir.ljust(25), 'structure created')


            # If VASP has not been run yet, we then create the job script for
            # the SLURM job scheduler.
            #if len(job_type)<4:

            if not isfile(struct_dir + outputfile):

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
                    os.system(job_command+' '+job_script_name)
                    os.chdir(twd)
                    print(struct_dir.ljust(25), 'job submitted')
                else:
                    print(struct_dir.ljust(25), 'job not yet run')
            # If VASP has already been run, we can write the resulting ionic
            # steps to ase trajectory files.
            else:
                if isfile(struct_dir + 'images.traj')==False or rebuild_traj_cache:

                    #what if image writting is interupred? taken care of by function 
                    images = convert_to_traj(struct_dir + outputfile, struct_dir + 'images.traj')	

                else:

                    size = os.path.getsize(struct_dir + 'images.traj')
                    #print(size)
                    if size > 4:
                        images = io.trajectory.Trajectory( struct_dir + 'images.traj', mode = 'r')
                    else:
                        images = []
                        print('%s was empty.'%(struct_dir + 'images.traj') )
                print(struct_dir.ljust(25), '  has {:4d} images'.format(len(images)))
                job_type_total_images += len(images)
        
        jtname = job_type[1] +'_' + job_type[2]
        print(jtname + ' Subtotal images: %i\n'% job_type_total_images)
        total_images += job_type_total_images
    print('Total images: %i'%total_images)


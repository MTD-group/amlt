




def vasp_job_maker(name_prefix, jobs, job_command, job_script_name, job_script_template, md_temperature_range = (100, 600), submit = False,
	random_structure_parameters = None, 
	known_structures=[], polymorphD3_parameters = None,
	first_structure = 'POSCAR.initial', magmom_filename = 'MAGMOMS.initial'):

	import numpy as np
	from ase import io

	from os.path import isfile
	import os
	from amlt import reasonable_random_structure_maker, polymorphD3, try_mkdir
	
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
				job_name = "%s_%s_%s_%i"%(name_prefix, job_type[1], job_type[2], structure_number)
				fid.write(job_script_template % (job_name, job_type[2])) 
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








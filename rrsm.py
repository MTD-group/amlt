

def try_mkdir(structure_direct):
	from os import mkdir
	try:
		mkdir(structure_direct)
	except:
		pass


from ase.data import covalent_radii as ase_covalent_radii
from numpy import array as np_array




def reasonable_random_structure_maker(elements, composition_generator,
										cut_off_radius = 5.0,
										fill_factor_max = 0.65, 
										fill_factor_min = 0.2,
										insert_attempt_max = 2000,
										max_build_failures = 10,
										element_radii = np_array(ase_covalent_radii),
										hard_radii    = 0.8*np_array(ase_covalent_radii),
										verbose = False):
	from ase import Atoms, Atom
	from ase.data import atomic_numbers, chemical_symbols
	from ase.geometry import find_mic
	z_numbers = [atomic_numbers[el] for el in elements ]
	
	from numpy import zeros, array, round, pi
	#argsort, floor,
	from numpy.random import rand

	print('Not tested with non-orthrombic systems!')

	def compute_nspecies(atomic_fractions, fill_factor = 1.0):
		atomic_volumes=[4/3.*pi*(element_radii[z])**3 for z in z_numbers]
		average_atomic_volume = array(atomic_fractions).dot(atomic_volumes)


		volume = (2*cut_off_radius)**3

		natoms = fill_factor * volume/average_atomic_volume

		nspecies = round(natoms *  atomic_fractions).astype(int )
		return nspecies


	def safe_insertion_test(atoms, new_species, position):
		if len(atoms)>0:
			safe = True
			D = atoms.get_positions() - position
			D_min, D_min_len = find_mic(D, cell =  atoms.get_cell())
			for existing_atom_index in range(len(atoms)):
				distance_cut = hard_radii[new_species] + \
						hard_radii[atoms.numbers[existing_atom_index]]
				if D_min_len[existing_atom_index] < distance_cut:
					safe = False
			return safe
		else:
			return True



	if verbose:
		print("Pure Compositions Possible:")
		for i in range(len(elements)):
			atomic_fractions = zeros(len(elements))
			atomic_fractions[i] = 1.0
			nspecies = compute_nspecies(atomic_fractions)
			print(elements[i].ljust(2), ':', int(nspecies.sum()))

		print('-------------------------------')



	composition_failed = True
	while composition_failed:

		fill_factor = (fill_factor_max - fill_factor_min) * rand() + fill_factor_min
		atomic_fractions =  composition_generator(elements)
		nspecies = compute_nspecies(atomic_fractions, fill_factor = fill_factor )
		
		last_prefix = 'Atom Count '+ ('(%i) '%sum(nspecies))
		
		print('Fill Factor'.ljust(len(last_prefix))+':', fill_factor)
		print('Composition'.ljust(len(last_prefix))+':', elements)
		print('Fractions'.ljust(len(last_prefix))  +':', atomic_fractions)
		print(last_prefix +':', nspecies)



		build_failed = True
		n_build_failures = 0
		while build_failed and n_build_failures < max_build_failures:
			atoms = Atoms(cell=[2*cut_off_radius, 2*cut_off_radius, 2*cut_off_radius], pbc = True)

			for species_index in range(len(nspecies)):
				for i in range(nspecies[species_index]):

					### multiple attempts per atom
					safe_insert = False
					attempt_no = 0

					while safe_insert==False and attempt_no < insert_attempt_max:
						position = atoms.get_cell().dot(rand(3))

						safe_insert = safe_insertion_test(atoms, z_numbers[species_index], position)
						attempt_no +=1
						#print(attempt_no)
					if safe_insert:
						atoms.append(Atom( elements[species_index], position) )

			#print(len(atoms), sum(nspecies))
			if len(atoms) == sum(nspecies):
				composition_failed = False
				build_failed = False
			else:
				n_build_failures += 1
				print('%i Build Failures'%n_build_failures)


	#print(structure_direct, 'Structure Successfully Made.')
	print('Structure Successfully Made.')
	
	return atoms












if __name__ == "__main__":



	elements = [ 'O', 'Si', 'Na', 'Ca']

	def generate_random_silicate(elements, no_gas = True):
		from numpy.random import rand
		atomic_fractions = rand(len(elements))
		atomic_fractions = atomic_fractions/atomic_fractions.sum()

		if no_gas:
			while atomic_fractions[0] > 0.70:
				atomic_fractions = rand(len(elements))
				atomic_fractions = atomic_fractions/atomic_fractions.sum()

		return atomic_fractions




	n_structures = 5
	for structure_number in range(n_structures):

		
		print('\n',structure_number)
		atoms = reasonable_random_structure_maker(elements, 
						fill_factor_max = 0.3, # 0.65 is about the max
						composition_generator = generate_random_silicate)

		print(atoms)

		if False:
			from ase import io
			from os.path import isfile
			try_mkdir(structure_direct)
			io.write(structure_direct+'POSCAR', atoms, format = 'vasp')


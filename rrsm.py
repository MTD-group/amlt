

def try_mkdir(structure_direct):
	from os import mkdir
	try:
		mkdir(structure_direct)
	except:
		pass


from ase.data import covalent_radii as ase_covalent_radii
from ase.data import chemical_symbols




default_element_radii = {}
default_hard_radii    = {}
for z in range(len(ase_covalent_radii)):
	default_element_radii[chemical_symbols[z]] = ase_covalent_radii[z]
	default_hard_radii[chemical_symbols[z]] = 0.9*ase_covalent_radii[z]

def reasonable_random_structure_maker(elements, composition_generator,
										cut_off_radius = 5.0,
										fill_factor_max = 0.55, 
										fill_factor_min = 0.2,
										insert_attempt_max = 2000,
										max_build_failures = 20,
										element_radii = default_element_radii,
										hard_radii    = default_hard_radii,
										magmom_generator = None,
										verbose = False):
	from ase import Atoms, Atom
	from ase.data import atomic_numbers, chemical_symbols
	from ase.geometry import find_mic
	#z_numbers = [atomic_numbers[el] for el in elements ]
	
	from numpy import zeros, array, round, pi
	#argsort, floor,
	from numpy.random import rand
	from random import shuffle

	print('Not tested with non-orthrombic systems!')

	def compute_nspecies(atomic_fractions, fill_factor = 1.0):
		atomic_volumes=[4/3.*pi*(element_radii[el])**3 for el in elements]
		average_atomic_volume = array(atomic_fractions).dot(atomic_volumes)


		volume = (2*cut_off_radius)**3

		natoms = fill_factor * volume/average_atomic_volume

		nspecies = round(natoms *  atomic_fractions).astype(int )
		return nspecies

	def compute_actual_fill_factor(atomic_fractions, nspecies):
		volume = (2*cut_off_radius)**3		

		atomic_volumes=[4/3.*pi*(element_radii[el])**3 for el in elements]
		
		total_atom_volume = array(atomic_volumes).dot(nspecies)

		
		return total_atom_volume/volume


	def safe_insertion_test(test_atoms, new_species, position):
		if len(test_atoms)>0:
			safe = True
			vec = test_atoms.get_positions() - position
			vec_min_image, vec_min_image_mag = find_mic(vec, cell =  test_atoms.get_cell())
			symbols = test_atoms.get_chemical_symbols()

			existing_atom_index = 0
			while existing_atom_index < len(test_atoms) and safe:
				distance_cut = hard_radii[new_species] + hard_radii[symbols[existing_atom_index]]

				if vec_min_image_mag[existing_atom_index] < distance_cut: 
					safe = False
				existing_atom_index += 1

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


		actual_fill_factor = compute_actual_fill_factor(atomic_fractions, nspecies)

		symbol_list = []
		for species_index in range(len(nspecies)):
			for i in range(nspecies[species_index]):
				symbol_list.append(elements[species_index])

		#print(symbol_list)

		output_positions = zeros((sum(nspecies),3))

		# random test order
		test_order = list(range(sum(nspecies)))
		shuffle(test_order)
		


		
		
		last_prefix = 'Atom Count '+ ('(%i) '%sum(nspecies))
		
		print('Fill Factor'.ljust(len(last_prefix))+': %.4f [%.4f]'% (actual_fill_factor,fill_factor))
		print('Composition'.ljust(len(last_prefix))+':', elements)
		print('Fractions'.ljust(len(last_prefix))  +':', atomic_fractions)
		print(last_prefix +':', nspecies)

		if True:
			test_symbols = ''
			for test_index in test_order:
				test_symbols = test_symbols + ' ' + symbol_list[test_index]
			print('Test Order'.ljust(len(last_prefix))  +':'+test_symbols)
                

		build_failed = True
		n_build_failures = 0
		
		while build_failed and n_build_failures < max_build_failures:

			# i use this to store only the safe atoms for tests, I rebuild so that element order is preserved later
			test_atoms = Atoms(cell=[2*cut_off_radius, 2*cut_off_radius, 2*cut_off_radius], pbc = True)


			safe_inserts = 0
			i = 0 
			attempt_no = 0
			while attempt_no < insert_attempt_max and i < sum(nspecies):
				test_index = test_order[i]
				### multiple attempts per atom
				safe_insert = False
				

				el = symbol_list[test_index]
				
			
				attempt_no = 0
				while safe_insert==False and attempt_no < insert_attempt_max:
					position = test_atoms.get_cell().dot(rand(3))

					safe_insert = safe_insertion_test(test_atoms, el, position)

					attempt_no +=1
					#print(attempt_no)
				if safe_insert:
					safe_inserts += 1
					test_atoms.append(Atom( el, position) )
					output_positions[test_index] = position

				i+=1

			#print(len(atoms), sum(nspecies))
			if safe_inserts == sum(nspecies):
				composition_failed = False
				build_failed = False

				atoms = Atoms(cell=[2*cut_off_radius, 2*cut_off_radius, 2*cut_off_radius], pbc = True)
				for i in range(sum(nspecies)):
					atoms.append(Atom( symbol_list[i], output_positions[i]) )

			else:
				n_build_failures += 1
				print('%i Build Failures'%n_build_failures)

	
	if callable(magmom_generator):
		magmoms = magmom_generator(atoms) #should return random spins for the atoms input
		atoms.set_initial_magnetic_moments(magmoms)
		print('Net Magnetic Moment:', magmoms.sum())
		
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
						fill_factor_max = 0.40, # 0.65 is about the max
						fill_factor_min = 0.2, #0.2 is default
						composition_generator = generate_random_silicate)

		#print(atoms)

		if False:
			from ase import io
			from os.path import isfile
			try_mkdir(structure_direct)
			io.write(structure_direct+'POSCAR', atoms, format = 'vasp')

		if True:
			from ase import io
			io.write('%i.POSCAR'%structure_number, atoms, format = 'vasp')


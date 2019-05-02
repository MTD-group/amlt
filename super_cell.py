#!/usr/bin/env python3




def compute_super_cell_needed_for_rcut(atoms, rcut):
	import numpy as np	
	cell = atoms.get_cell()
	recip_cell = atoms.get_reciprocal_cell()	
	cells_needed = [1,1,1]
	for i in range(3):
		scaled_kvec = [0,0,0]
		scaled_kvec[i] = 1.0
		kvec = recip_cell.dot(scaled_kvec)
		spacing = 1.0/(np.sqrt(kvec.dot(kvec)))
		cells_needed[i] = int(   np.ceil( 2.0*rcut/ spacing )  )

	return cells_needed

def super_cell(atoms, cells_needed, use_initial_magnetic_moments = False):
	from ase import io, Atoms, Atom
	from ase.calculators.singlepoint import SinglePointCalculator
	import numpy as np

	cell     = atoms.get_cell()
	elements = atoms.get_chemical_symbols()
	forces   = atoms.get_forces()
	energy   = atoms.get_potential_energy()


	cell_out = []
	for dim in range(3):
		cell_out.append(cells_needed[dim]*cell[dim])
	atoms_out = Atoms(pbc= True, cell = cell_out)
	forces_out  = []
	energy_out  = 0.0

	# do a test for magmoms
	try:	
		magmoms  = atoms.get_magnetic_moments()
		magmoms_out = []
	except:
		if use_initial_magnetic_moments:
			magmoms = atoms.get_initial_magnetic_moments()
			magmoms_out = []
		else:	
			magmoms = len(atoms)*[None]
			magmoms_out = None

	

	for i in range(cells_needed[0]):
		for j in range(cells_needed[1]):
			for k in range(cells_needed[2]):
				energy_out += energy
				shift = cell.dot([i,j,k])
				for atom_index in range(len(atoms)):
					atoms_out.append(Atom(	elements[atom_index],
											position = atoms.positions[atom_index] + shift,
											magmom = magmoms[atom_index],
											charge = None))
					forces_out.append(forces[atom_index])

					if magmoms_out != None:
						magmoms_out.append(magmoms[atom_index])
	
	calc = SinglePointCalculator(atoms_out,
			                     energy=energy_out,
			                     forces=forces_out,
			                     magmoms =magmoms_out)
	atoms_out.set_calculator(calc)

	return atoms_out



def super_cell_if_needed(atoms, rcut, verbose = True, use_initial_magnetic_moments = False):
	cells_needed = compute_super_cell_needed_for_rcut(atoms, rcut)
	if sum(cells_needed) == 3:
		return atoms
	else:
		if verbose: 
			print ("Stucture expanded to %ix%ix%i" %tuple(cells_needed) + " super cell to fit rcut = %f"%rcut)
		return super_cell(atoms, cells_needed, use_initial_magnetic_moments = use_initial_magnetic_moments )





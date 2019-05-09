
from .super_cell import compute_super_cell_needed_for_rcut
from ase import build
import numpy as np

# our version of the super_cell should only be used for padding outputs up to the rcut for AMP
def random_super_cell(atoms, rcut):
	"""
	Args:
		atoms (Atoms):  input structure
		rcut (float): maximum distance for considering pairs of atoms
	Returns:
		atoms_out (Atoms): new super cell
	"""
	
	cell_ranges = compute_super_cell_needed_for_rcut(atoms, rcut = rcut)
	
	new_cells = []
	for i in range(3):
		line = [0,0,0]
		line[i] = np.random.random_integers(1,cell_ranges[i])
		new_cells.append(line)
	
	atoms_out = build.cut(atoms, a = new_cells[0], b = new_cells[1], c = new_cells[2] )
	
	return(atoms_out)
	
	
def random_deletion(atoms,  deletion_chance):
	""" Add deletions to atoms object in place.
	
	Args:
		atoms (Atoms): ASE atoms object to modify
		deletion_chance (Float): fraction of atoms to remove on average
	"""
	'''Deletion chance is fraction'''
	mask = np.random.rand(len(atoms))
	mask = np.where(mask < deletion_chance, True, False)
	del atoms[mask]
	
	
def random_distortion(atoms, atom_distortion, lattice_distortion, volume_change_max = 0.05):
	""" Add Gaussian noise to atoms's positions and lattice in place.
	Args:
		atoms (Atoms): ASE atoms object to modify
		atom_distortion (Float):  variance of Gaussian to add
			to atom positions in Angstroms
		lattice distortion (Float): variance of Gaussian to
			multiply cell lengths and angles by in percent
						
	"""
	scaled_positions = atoms.get_scaled_positions()
	vol0 = atoms.get_volume()
	old_cell =  atoms.get_cell()
	old_lengths = atoms.get_cell_lengths_and_angles()[0:3]
	volume_change = 2.0*volume_change_max # just to start the loop	
	next_cell = np.zeros((3,3))

	while abs(volume_change)>= abs(volume_change_max):
		for dir_index in range(3):
			l = old_lengths[dir_index]
			delta = np.random.randn(3)
			delta = delta/np.sqrt(delta.dot(delta))
			next_cell[dir_index] = lattice_distortion*l*delta +  old_cell[dir_index]

		atoms.set_cell( next_cell )
		volume_change = (atoms.get_volume()-vol0)/vol0		

	atoms.set_scaled_positions(scaled_positions) # this prevents missalignment of atoms with a radically shifted cell
	atoms.wrap()

	## this is the old unlreliable method
	## new_cell = (1.0+ lattice_distortion*np.random.randn(6)) * atoms.get_cell_lengths_and_angles()
	## print( new_cell)
	
	
	# small atom distortions
	atoms.set_positions ( atoms.positions + atom_distortion* np.random.randn(len(atoms),3) )

def random_magnetic_moment_flips(atoms, flip_chance):
	mask = np.random.rand(len(atoms))
	mask = np.where(mask < flip_chance, True, False)

	magmoms = atoms.get_initial_magnetic_moments()
	magmoms[mask] = -magmoms[mask]
	
	atoms.set_initial_magnetic_moments(magmoms)
	
def polymorphD3(atoms, atom_distortion=0.2, lattice_distortion=0.10, deletion_chance=0.05, rcut=6.5, volume_change_max = 0.05, flip_chance = 0.10): 
	""" Creates a perturbed version of the input structure.
	
	Creates a random supercell, removes random atoms, and randomly 
	distorts atom positions and lattice.
	
	Args:
		atoms (Atoms): ASE atoms object
		atom_distortion (Float):  variance of Gaussian to add
			to atom positions in Angstroms
		lattice distortion (Float): variance of Gaussian to
			multiply cell lengths and angles by in percent
		deletion_chance (Float): fraction of atoms to remove on average
		rcut (float): maximum distance for considering pairs of atoms
	Returns:
		atoms_out (Atoms): new perturbed structure
	"""
	atoms_out = random_super_cell(atoms, rcut)
	
	random_deletion(atoms_out, deletion_chance)
	
	random_distortion(atoms_out, atom_distortion, lattice_distortion, volume_change_max)
 
	random_magnetic_moment_flips(atoms, flip_chance)

	return atoms_out

if __name__ == "__main__":

	from ase import io
	name = 'TlInS2_mp-632539_primitive.cif'
	unit_cell = io.read('test_structures/'+'TlInS2_mp-632539_primitive.cif')
	
	atoms = polymorphD3(unit_cell,  rcut = 5.0 )
	
	name_out = name.replace('.cif', '.poscar' )
	io.write('test_results/'+name_out, atoms)

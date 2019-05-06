
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
	
	
def random_distortion(atoms, atom_distortion, lattice_distortion):
	""" Add Gaussian noise to atoms's positions and lattice in place.
	Args:
		atoms (Atoms): ASE atoms object to modify
		atom_distortion (Float):  variance of Gaussian to add
			to atom positions in Angstroms
		lattice distortion (Float): variance of Gaussian to
			multiply cell lengths and angles by in percent
						
	"""
	
	new_cell = (1.0+ lattice_distortion*np.random.randn(6)) * atoms.get_cell_lengths_and_angles()
	atoms.set_cell( new_cell )
	atoms.set_positions ( atoms.positions + atom_distortion* np.random.randn(len(atoms),3) )

	
def polymorphD3(atoms, atom_distortion=0.2, lattice_distortion=0.10, deletion_chance=0.05, rcut=6.5): 
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
	
	random_distortion(atoms_out, atom_distortion, lattice_distortion)
 
	return atoms_out

if __name__ == "__main__":

	from ase import io
	name = 'TlInS2_mp-632539_primitive.cif'
	unit_cell = io.read('test_structures/'+'TlInS2_mp-632539_primitive.cif')
	
	atoms = polymorphD3(unit_cell,  rcut = 5.0 )
	
	name_out = name.replace('.cif', '.poscar' )
	io.write('test_results/'+name_out, atoms)

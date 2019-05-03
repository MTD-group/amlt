
from super_cell import compute_super_cell_needed_for_rcut
from ase import build
import numpy as np

# our version of the super_cell should only be used for padding outputs up to the rcut for AMP
def random_super_cell(atoms, rcut):

	cell_ranges = compute_super_cell_needed_for_rcut(atoms, rcut = rcut)
	
	#print(cell_ranges)
	new_cells = []
	for i in range(3):
		line = [0,0,0]
		line[i] = np.random.random_integers(1,cell_ranges[i])
		new_cells.append(line)
	
	#print(new_cells)
	atoms_out = build.cut(atoms, a = new_cells[0], b = new_cells[1], c = new_cells[2] )
	
	return(atoms_out)
	
def random_deletion(atoms,  deletion_chance):
	'''Deletion chance is fraction'''
	mask = np.random.rand(len(atoms))
	mask = np.where(mask < deletion_chance, True, False)
	#print(mask)
	del atoms[mask]
	
def random_distortion(atoms, atom_distortion, lattice_distortion):
	'''atom_distortion is in Ang and lattice distortion is in percent'''
	
	#print(atoms.get_cell_lengths_and_angles())
	new_cell = (1.0+ lattice_distortion*np.random.randn(6)) * atoms.get_cell_lengths_and_angles()
	#print(new_cell)
	atoms.set_cell( new_cell )
	#print(atoms.get_positions()[0])
	atoms.set_positions ( atoms.positions + 0.1* np.random.randn(len(atoms),3) )
	#print(atoms.get_positions()[0])

def polymorphD3(atoms, atom_distortion = 0.2, lattice_distortion = 0.10, deletion_chance=0.05, rcut = 6.5): 

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

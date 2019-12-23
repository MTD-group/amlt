
from super_cell import compute_super_cell_needed_for_rcut
from ase import build
import numpy as np



    
    
class PolymorphD3(object):
    
    def __init__(self, 
                 atoms,
		 elements = None, 
                 atom_distortion=0.2, 
                 lattice_distortion=0.10,
		 shrink_bias = 0.25, 
                 deletion_chance=0.05, 
                 rcut=6.5, 
                 volume_change_max = 0.05, 
                 flip_chance = 0.10,
		 swap_chance = 0.05,
                 random_seed = None): 
    
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
	    swap_chance (Float): fraction of atoms to switch element type
            rcut (Float): maximum distance for considering pairs of atoms
            volume_change_max (Float): Relative amount volume allowed to change
            flip_chance (Float): If None, no flips; otherwise chance of 
                flipping a magnetic moment's direction
            random_seed (Integer): Set for reproducible Numpy randomness
        """
        self.atoms = atoms
        if isinstance(random_seed, int):
            np.random.seed(random_seed)
        self.atoms_out = self.random_super_cell(self.atoms, rcut)
        
        self.random_deletion(self.atoms_out, deletion_chance)
        
        self.random_distortion(self.atoms_out, 
                          atom_distortion, 
                          lattice_distortion, 
                          volume_change_max,
			  shrink_bias)
			  
        if elements == None:
            el = atoms.get_chemical_symbols()
            elements = list(set(el))
	    
        self.random_swaps(self.atoms_out, elements, swap_chance)
        
        if flip_chance:
            self.random_magnetic_moment_flips(self.atoms_out, flip_chance)
        
    
    def random_super_cell(self, atoms, rcut):
        """
        Args:
            atoms (Atoms):  input structure
            rcut (float): maximum distance for considering pairs of atoms
        Returns:
            atoms_out (Atoms): new super cell
            
        Our version of the super_cell should only be used for padding outputs 
        up to the rcut for AMP.
        """
        
        cell_ranges = compute_super_cell_needed_for_rcut(atoms, 
                                                         rcut = rcut)
        
        new_cells = []
        for i in range(3):
            line = [0,0,0]
            line[i] = np.random.randint(1,cell_ranges[i])
            new_cells.append(line)
        
        atoms_out = build.cut(atoms, 
                              a = new_cells[0], 
                              b = new_cells[1], 
                              c = new_cells[2] )
        return(atoms_out)
    
    def random_deletion(self, atoms,  deletion_chance):
        """ Add deletions to atoms object in place.
        
        Args:
            atoms (Atoms): ASE atoms object to modify
            deletion_chance (Float): fraction of atoms to remove on average
        """
        '''Deletion chance is fraction'''
        mask = np.random.rand(len(atoms))
        mask = np.where(mask < deletion_chance, True, False)
        del atoms[mask]
	
    def random_swaps(self, atoms, elements, swap_chance):
        mask = np.random.rand(len(atoms))
        mask = np.where(mask < swap_chance, True, False)
        el = atoms.get_chemical_symbols()
	
        if len(elements)>1:
            for ia in range(len(atoms)):
                if mask[ia]:
                    ie = np.random.randint(0,len(elements))
                    el[ia] = elements[ie]
            atoms.set_chemical_symbols(el)
    
    def random_distortion(self, 
                          atoms, 
                          atom_distortion, 
                          lattice_distortion, 
                          volume_change_max = 0.05,
			  shrink_bias = 0.25):
        """ Add Gaussian noise to atoms's positions and lattice in place.
        Args:
            atoms (Atoms): ASE atoms object to modify
            atom_distortion (Float):  variance of Gaussian to add
                to atom positions in Angstroms
            lattice distortion (Float): variance of Gaussian to
                multiply cell lengths and angles by in percent
            volume_change_max (Float): Relative amount volume allowed to change
                            
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
                delta = np.random.randn(3)-0.5
                delta = delta/np.sqrt(delta.dot(delta))
                next_cell[dir_index] = (lattice_distortion*l*delta + old_cell[dir_index])
    
            atoms.set_cell( next_cell )
            volume_change = (atoms.get_volume()-vol0)/vol0        
    
        dilation = lattice_distortion*(np.random.rand()-shrink_bias) + 1.0
        atoms.set_cell(dilation*atoms.get_cell())
	
	# this prevents missalignment of atoms with a radically shifted cell
        atoms.set_scaled_positions(scaled_positions) 
        atoms.wrap()
        
        # small atom distortions
        atoms.set_positions(atoms.positions + 
                            atom_distortion * 
                            np.random.randn(len(atoms),3) )
    
    def random_magnetic_moment_flips(self, atoms, flip_chance):
        """Randomly mutates direction of spins."""
        mask = np.random.rand(len(atoms))
        mask = np.where(mask < flip_chance, True, False)
    
        magmoms = atoms.get_initial_magnetic_moments()
        magmoms[mask] = -magmoms[mask]
        
        atoms.set_initial_magnetic_moments(magmoms)
    
        

if __name__ == "__main__":

    from ase import io
    name = 'TlInS2_mp-632539_primitive.cif'
    unit_cell = io.read('test_structures/'+'TlInS2_mp-632539_primitive.cif')
    
    Poly = PolymorphD3(unit_cell,  rcut = 5.0, flip_chance=0.10, random_seed=0)
    atoms = Poly.atoms_out
    
    name_out = name.replace('.cif', '.poscar' )
    io.write('test_results/'+name_out, atoms)

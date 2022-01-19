
from .super_cell import compute_super_cell_needed_for_rcut
from .random_sheer import random_sheer_matrix

from ase import Atoms
import numpy as np


rng = np.random.default_rng() 
    


def random_cell_strain(
                        cell,
                        volume_strain=0.10,
                        von_mises_strain=0.20,
                        shrink_bias = 0.25,
                        rng=rng):
    """ .
    Args:
        atoms (Atoms): ASE atoms object to modify
        atom_distortion (Float):  variance of Gaussian to add
            to atom positions in Angstroms
        lattice distortion (Float): variance of Gaussian to
            multiply cell lengths and angles by in percent
        volume_change_max (Float): Relative amount volume allowed to change
                        
    """

    vms = von_mises_strain*rng.random()
    M = random_sheer_matrix(vms,rng=rng)

    dilation =volume_strain*(rng.random()-shrink_bias) + 1.0
    M = dilation*M
    cell_new = M @ cell
    return cell_new




def random_super_cell(cell, rcut, min_cells = 2, rng=rng):
    """
    Args:
        atoms (Atoms):  input structure
        rcut (float): maximum distance for considering pairs of atoms
        min_cells (tuple/list of ints): the minimum number of cells ( = N1 x N2 x N3) to include in the super_cell
    Returns:
        atoms_out (Atoms): new super cell
        
    Our version of the super_cell should only be used for padding outputs 
    up to the rcut for AMP.
    """
    
    cell_max = compute_super_cell_needed_for_rcut(cell=cell, rcut = rcut)
    rint = rng.integers #(low, high=None, size=None )# high is 1+      

    cells = np.ones(3,dtype=int)
    for i in range(3):
        cells[i] = rint(low=1, high=cell_max[i]+1, size=1) 
    
    # now we expand the duplication till we hit the min_cells
    while cells[0]*cells[1]*cells[2] < min_cells:
        randdim = rint(low=1,high=3, size=1)
        cells[randdim]+=1

    return(cells)



def random_deletion( atoms, deletion_chance=0.05, rng=rng):
    """ Add deletions to atoms object in place.
    
    Args:
        atoms (Atoms): ASE atoms object to modify
        deletion_chance (Float): fraction of atoms to remove on average
    """
    '''Deletion chance is fraction'''
    mask = rng.random(len(atoms))
    mask = np.where(mask < deletion_chance, True, False)
    del atoms[mask]
    




def random_element_swaps(atoms, elements=None, swap_chance=0.05, rng=rng):

    el = atoms.get_chemical_symbols()
    if elements is None:
        elements = list(set(el))
    if len(elements)>1:
        mask = np.where(rng.random(len(atoms)) < swap_chance, True, False)
        for ia in range(len(atoms)):
            if mask[ia]:
                ie = rng.integers(0,len(elements))
                el[ia] = elements[ie]
        atoms.set_chemical_symbols(el)



def random_magnetic_moment_flips(atoms, flip_chance=0.10, rng=rng):
    """Randomly mutates direction of spins."""
    mask = np.where(rng.random(len(atoms)) < flip_chance, True, False)

    magmoms = atoms.get_initial_magnetic_moments()
    magmoms[mask] = -magmoms[mask]
    
    atoms.set_initial_magnetic_moments(magmoms)



def polymorphate(
                 atoms,
                 elements = None, 
                 rcut=6.0,
                 atom_displacement=0.2, 
                 volume_strain=0.10,
                 von_mises_strain = 0.20, 
                 shrink_bias = 0.25, 
                 deletion_chance=0.05, 
                 flip_chance = 0.10,
                 swap_chance = 0.05,
                 min_cells = 2,
                 rng = rng):


        """ Creates a perturbed version of the input structure.
        
        Creates a random supercell, removes random atoms, and randomly 
        distorts atom positions and lattice.
        
        Args:
            atoms (Atoms): ASE atoms object
            atom_distortion (Float):  variance of Gaussian to add
                to atom positions in Angstroms
            volume_strain (Float): range of random strain to
                multiply cell lengths and angles by in percent
            deletion_chance (Float): fraction of atoms to remove on average
        swap_chance (Float): fraction of atoms to switch element type
            rcut (Float): maximum distance for considering pairs of atoms
            volume_change_max (Float): Relative amount volume allowed to change
            flip_chance (Float): If None, no flips; otherwise chance of 
                flipping a magnetic moment's direction
            random_seed (Integer): Set for reproducible Numpy randomness
        """


        atoms_out = atoms.copy()

        ##### Strain the cell
        new_cell = random_cell_strain(
                            cell=atoms_out.cell,
                            volume_strain=volume_strain,
                            von_mises_strain = von_mises_strain,
                            shrink_bias=shrink_bias,
                            rng=rng)
        atoms_out.set_cell(new_cell, scale_atoms=True)

        #### super cell
        cells = random_super_cell(atoms_out.cell, rcut,  min_cells=min_cells, rng=rng)
        atoms_out = atoms_out.repeat(cells)
        
        
        #### deletions
        random_deletion(atoms_out, deletion_chance=deletion_chance, rng=rng)
        
        ### elements swaps
        random_element_swaps(atoms_out, elements=elements, swap_chance=swap_chance, rng=rng)
        
        
        ### displace atoms
        atoms_out.rattle(atom_displacement,rng=rng)
        
        ### magmoms
        if np.linalg.norm(atoms.get_initial_magnetic_moments()) > 1e-5:
            self.random_magnetic_moment_flips(atoms_out, flip_chance=flip_chance, rng=rng)
        
        return atoms_out


class Polymorpher2(object):
    
    def __init__(self, 
                 atoms,
                 weights=None,
                 elements = None, 
                 rcut=6.0,
                 atom_displacement=0.2, 
                 volume_strain=0.10,
                 von_mises_strain = 0.20, 
                 shrink_bias = 0.25, 
                 deletion_chance=0.05, 
                 flip_chance = 0.10,
                 swap_chance = 0.05,
                 min_cells = 2,
                 rng = rng): 
    
        """ Creates a perturbed version of the input structure.
        
        Creates a random supercell, removes random atoms, and randomly 
        distorts atom positions and lattice.
        
        Args:
            atoms (Atoms): ASE atoms object
            atom_distortion (Float):  variance of Gaussian to add
                to atom positions in Angstroms
            volume_strain (Float): range of random strain to
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
        if weights is None:
            self.weights=np.ones(len(atoms))/len(atoms)

        self.elements=elements
        self.rcut=rcut
        self.atom_displacement=atom_displacement
        self.volume_strain=volume_strain
        self.von_mises_strain=von_mises_strain
        self.shrink_bias=shrink_bias
        self.deletion_chance=deletion_chance
        self.flip_chance=flip_chance
        self.swap_chance=swap_chance
        self.min_cells=min_cells
        self.rng = rng
        
    def __call__(self, return_index = False):
        if type(self.atoms) is Atoms:
            atoms_picked = self.atoms
            index = 0
        else:
            index = self.rng.choice(len(self.atoms), p=self.weights/sum(self.weights))
            atoms_picked = self.atoms[index]

        atoms_out = polymorphate(
                atoms_picked,
                elements =self.elements, 
                rcut=self.rcut,
                atom_displacement=self.atom_displacement, 
                volume_strain=self.volume_strain,
                von_mises_strain =self.von_mises_strain, 
                shrink_bias = self.shrink_bias, 
                deletion_chance=self.deletion_chance, 
                flip_chance = self.flip_chance,
                swap_chance =self.swap_chance,
                min_cells = self.min_cells,
                rng=self.rng)

        if return_index:
            return atoms_out, index
        else:
            return atoms_out
    

    


    

    
        



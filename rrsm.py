

from ase import Atoms, Atom
from ase.cell import Cell
from ase.data import atomic_numbers, chemical_symbols
from ase.geometry import find_mic
import numpy as np 

rng = np.random.default_rng()


from ase.data import covalent_radii as ase_covalent_radii
from ase.data import chemical_symbols

default_hard_radius_ratio = 0.85




def get_hard_radii_by_ratio(element_radii,  hard_radius_ratio = default_hard_radius_ratio):
    hard_radii = {}
    for z in range(1,len(element_radii)):
        hard_radii[chemical_symbols[z]] = hard_radius_ratio*element_radii[chemical_symbols[z]]
    return hard_radii


default_element_radii = {}
for z in range(len(ase_covalent_radii)):
    default_element_radii[chemical_symbols[z]] = ase_covalent_radii[z]
default_hard_radii    = get_hard_radii_by_ratio(default_element_radii)


def compute_element_volumes(elements, element_radii):
    return [4/3.*np.pi*(element_radii[el])**3 for el in elements]
    
    
def compute_nspecies(atomic_fractions, element_volumes, cell_volume, fill_factor = 0.5 ):
    average_atomic_volume = np.array(atomic_fractions).dot(element_volumes)
    natoms = fill_factor * cell_volume/average_atomic_volume
    nspecies = np.ceil(natoms *  atomic_fractions).astype(int )
    return nspecies


def compute_actual_fill_factor(nspecies, element_volumes, cell_volume):
    total_atom_volume = np.array(element_volumes).dot(nspecies)
    return total_atom_volume/cell_volume


def generate_random_symbol_list(
                                composition_generator,  
                                fill_factor_min, fill_factor_max, 
                                elements, element_volumes, 
                                cell_volume,  
                                rng=rng,
                                return_details = True):
                                
    fill_factor = (fill_factor_max - fill_factor_min) * rng.random() + fill_factor_min 
    atomic_fractions =  composition_generator(rng=rng)
    nspecies = compute_nspecies(atomic_fractions, element_volumes, cell_volume, fill_factor = fill_factor )

    symbol_list = []
    for species_index in range(len(nspecies)):
        for i in range(nspecies[species_index]):
            symbol_list.append(elements[species_index])
    
    rng.shuffle(symbol_list)
    
    if return_details:
        return symbol_list, fill_factor, atomic_fractions, nspecies
    else:
        return symbol_list


def safe_insertion_test(test_atoms, new_species, position, hard_radii):
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


def reasonable_random_structure_maker(elements, 
                                        composition_generator=None,
                                        cell =12.0,
                                        fill_factor_max = 0.50, 
                                        fill_factor_min = 0.05,
                                        insert_attempt_max = 1000,
                                        max_build_failures = 20,
                                        element_radii = default_element_radii,
                                        hard_radii    = None,
                                        hard_radius_ratio = default_hard_radius_ratio,
                                        magmom_generator = None,
                                        verbose = False,
                                        max_atoms = 200,
                                        rng=rng,
                                        ):
    
    if hard_radii is None:
        hard_radii =  get_hard_radii_by_ratio(element_radii, default_hard_radius_ratio)


    if composition_generator is None:
        ntypes=len(elements)
        def composition_generator(rng):
            comp =np.absolute( rng.normal(size=ntypes))
            return comp/comp.sum()

    if not hasattr(cell, '__iter__'):
        cell = 3*[float(cell)]
        
    atoms=Atoms(cell=cell)# this temporary atoms object lets us be more flxible in the cell format
    cell_volume = atoms.cell.volume
    simcell = atoms.get_cell()

    element_volumes = compute_element_volumes(elements, element_radii)
    

    if verbose:
        print("Pure Compositions Possible at %.3f fill_factor:"%fill_factor_max)
        for i in range(len(elements)):
            atomic_fractions = np.zeros(len(elements))
            atomic_fractions[i] = 1.0
            nspecies = compute_nspecies(atomic_fractions, element_volumes, cell_volume, fill_factor = fill_factor_max )
            print(elements[i].ljust(2), ':', int(nspecies.sum()))
        print('----------------------------------')


    ############ now we attempt to build structures
    composition_failed = True
    while composition_failed:
        natoms = max_atoms+1
        while natoms > max_atoms:
            symbol_list, fill_factor, atomic_fractions, nspecies = generate_random_symbol_list(
                                            composition_generator,  
                                            fill_factor_min, 
                                            fill_factor_max,
                                            elements,
                                            element_volumes, 
                                            cell_volume,  
                                            rng=rng,
                                            return_details=True)
            natoms = len(symbol_list)
            
        if verbose:
            actual_fill_factor = compute_actual_fill_factor(nspecies, element_volumes, cell_volume)
            last_prefix = 'Atom Count '+ ('(%i) '%sum(nspecies))
            print('Fill Factor'.ljust(len(last_prefix))+': %.4f [%.4f]'% (actual_fill_factor,fill_factor))
            print('Composition'.ljust(len(last_prefix))+':', elements)
            print('Fractions'.ljust(len(last_prefix))  +':', atomic_fractions)
            print(last_prefix +':', nspecies)

            if True:
                symbols_order_string = ''
                for sym in symbol_list:
                    symbols_order_string += ' '+sym
                print('Test Order'.ljust(len(last_prefix))  +':'+symbols_order_string)
                

        build_failed = True
        n_build_failures = 0
        
        while build_failed and n_build_failures < max_build_failures:
            # if all the atoms are inserted into this structure, then it's complete
            atoms = Atoms(cell=simcell, pbc = True)

            safe_inserts = 0
            atom_index = 0 
            attempt_no = 0
            while attempt_no < insert_attempt_max and atom_index < len(symbol_list):
                safe_insert = False
                new_species = symbol_list[atom_index]

                attempt_no = 0
                while safe_insert==False and attempt_no < insert_attempt_max:
                    position = rng.random(3) @ atoms.cell 
                    safe_insert = safe_insertion_test(atoms, new_species, position, hard_radii)
                    attempt_no +=1

                if safe_insert:
                    safe_inserts += 1
                    atoms.append(Atom( new_species, position) )

                atom_index+=1

            if safe_inserts == sum(nspecies):
                composition_failed = False
                build_failed = False
                
            else:
                n_build_failures += 1
                if verbose: 
                    print('%i Build Failures'%n_build_failures)

    
    if callable(magmom_generator):
        magmoms = magmom_generator(atoms,rng=rng) #should return random spins for the atoms input
        atoms.set_initial_magnetic_moments(magmoms)
        print('Net Magnetic Moment:', magmoms.sum())
        
    if verbose:
        print('Structure Successfully Made.')
    
    return atoms





class RandomStructure:

    def __init__(self,
            elements, 
            composition_generator=None,
            cell =12.0,
            fill_factor_max = 0.50, 
            fill_factor_min = 0.05,
            insert_attempt_max = 1000,
            max_build_failures = 20,
            element_radii = default_element_radii,
            hard_radii    = None,
            hard_radius_ratio = default_hard_radius_ratio,
            magmom_generator = None,
            verbose = False,
            max_atoms = 200,
            rng=rng):
        
        
        self.elements              = elements
        self.composition_generator = composition_generator
        self.cell                  = cell
        self.fill_factor_max       = fill_factor_max
        self.fill_factor_min       = fill_factor_min
        self.insert_attempt_max    = insert_attempt_max
        self.max_build_failures    = max_build_failures
        self.element_radii         = element_radii
        self.hard_radii            = hard_radii
        self.hard_radius_ratio     = hard_radius_ratio
        self.magmom_generator      = magmom_generator
        self.verbose               = verbose
        self.rng                   = rng
        self.max_atoms             = max_atoms
    
    def __call__(self):
        atoms = reasonable_random_structure_maker(
                                        elements=self.elements, 
                                        composition_generator = self.composition_generator,
                                        cell = self.cell,
                                        fill_factor_max = self.fill_factor_max, 
                                        fill_factor_min = self.fill_factor_min,
                                        insert_attempt_max = self.insert_attempt_max,
                                        max_build_failures = self.max_build_failures,
                                        element_radii = self.element_radii,
                                        hard_radii    = self.hard_radii,
                                        hard_radius_ratio = self.hard_radius_ratio,
                                        magmom_generator = self.magmom_generator,
                                        verbose = self.verbose,
                                        max_atoms = self.max_atoms,
                                        rng=self.rng,
                                        )

        return atoms





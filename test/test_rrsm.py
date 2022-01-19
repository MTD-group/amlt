

from amlt.rrsm import RandomStructure
from ase import io
import numpy as np 

rng = np.random.default_rng(8478381835)


elements = [ 'O', 'Si', 'Na', 'Ca']

def generate_random_silicate(rng, no_gas = True):
    atomic_fractions = np.absolute(rng.normal(size=len(elements)))
    atomic_fractions = atomic_fractions/atomic_fractions.sum()
    if no_gas:
        while atomic_fractions[0] > 0.70:
            atomic_fractions = np.absolute(rng.normal(size=len(elements)))
            atomic_fractions = atomic_fractions/atomic_fractions.sum()
    return atomic_fractions



randstruct = RandomStructure(
                    elements = elements, 
                    cell = 12.0,
                    fill_factor_max = 0.40, 
                    fill_factor_min = 0.10, 
                    verbose = True,
                    composition_generator = generate_random_silicate,
                    rng=rng)


###### cube test###
n_structures = 5
traj = io.Trajectory('random_structures.traj', 'w')
for structure_number in range(n_structures):
    print()
    print('->',structure_number)
    atoms = randstruct()

    traj.write(atoms)
traj.close()



#### non-orthorhombic test ####
nonboxcell = [[0,12.0,12],[12,0,12],[12,12,0]]
randstruct.cell= nonboxcell

n_structures = 5
traj = io.Trajectory('random_structures_non-orthorhombic.traj', 'w')
for structure_number in range(n_structures):
    print()
    print('->',structure_number)
    atoms = randstruct()

    traj.write(atoms)
traj.close()


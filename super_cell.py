#!/usr/bin/env python3
import numpy as np  



def compute_super_cell_needed_for_rcut(cell, rcut):
    rcell = cell.reciprocal()  
    cells_needed = np.ones(3)
    for i in range(3):
        spacing = 1.0/np.linalg.norm( rcell[i])
        cells_needed[i] = int(   np.ceil( 2.0*rcut/ spacing )  )
    return cells_needed


def super_cell(atoms, cells_needed, use_initial_magnetic_moments = False):
    from ase import io, Atoms, Atom
    from ase.calculators.singlepoint import SinglePointCalculator

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
    has_magmoms = False
    calc = atoms.get_calculator()
    if calc != None:
        if 'magmom' in calc.results:    
            has_magmoms = True
    
    if has_magmoms:
        magmoms  = atoms.get_magnetic_moments()
        magmoms_out = [] # will be filled below
    else:
        magmoms = len(atoms)*[None]
        magmoms_out = None
    
    # allows overriding info
    if use_initial_magnetic_moments:
        magmoms = atoms.get_initial_magnetic_moments()
        magmoms_out = []
        
            

    
    dup_index = np.zeros(3)
    for i in range(cells_needed[0]):
        dup_index[0] = i
        for j in range(cells_needed[1]):
            dup_index[1] = j
            for k in range(cells_needed[2]):
                dup_index[2] = k
                
                energy_out += energy
                shift = dup_index.dot(cell)
                for atom_index in range(len(atoms)):
                    atoms_out.append(Atom(    elements[atom_index],
                                            position = atoms.positions[atom_index] + shift,
                                            magmom = magmoms[atom_index],
                                            charge = None))
                    forces_out.append(forces[atom_index])

                    if magmoms_out != None:
                        magmoms_out.append(magmoms[atom_index])
    # calc.results dict can edited directly
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





if __name__ == "__main__":

    from ase import io
    
    traj = io.Trajectory('test_structures/AFM_Cr.traj','r')
    Cr_atoms = traj[0]
    traj.close()
    
    ZrO2_atoms = io.read('test_structures/ZrO2.OUTCAR')
    
    #tlins_atoms = io.read('test_structures/TlInS2_mp-632539_primitive.cif')

    tests = [Cr_atoms, ZrO2_atoms]

    for atoms in tests:
        super_cell(atoms, [1,1,2])
    
        #calc = atoms.get_calculator()
        #print(calc)
        #if calc != None:
        #    if 'magmom' in calc.results:
        #        atoms.get_magnetic_moments()

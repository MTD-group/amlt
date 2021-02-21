import pytest
from ase import Atoms
from amlt import contour_exploration
import numpy as np
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms


def test_kappa():
    
    radius = 2.6 
    atoms = Atoms(pair, positions=[[0, 0, 0],[0,radius,0]])
    atoms.set_constraint(FixAtoms(indices=[0]))
    atoms.center(vacuum=10)
    atoms.calc = EMT()

    atoms.set_calculator(EMT())
    initial_energy = atoms.get_potential_energy()



    print("Target Curvature {: .4f} eV/atoms".format( (1/radius)))

    name = 'test_kappa1'
    traj_name = name+'.traj'
    log_name = name+'.log'


    dyn = contour_exploration(atoms,
                    maxstep = 1.0,
                    parallel_drift = 0.0,
                    remove_translation  = False,
                    force_parallel_step_scale = None,
                    energy_target = initial_energy,
                    angle_limit = 20,
                    use_tangent_curvature= False,
                    #trajectory = traj_name,
                    #logfile = log_name,
                    ) 

    
    for i in range(5):
        dyn.run(10)
        energy_error = (atoms.get_potential_energy()-initial_energy)/len(atoms)
        print('Potentiostat Error {: .4f} eV/atom'.format( energy_error))
        assert 0 == pytest.approx(energy_error, abs=0.01)
        

if __name__ == '__main__':
    test_potentiostat()

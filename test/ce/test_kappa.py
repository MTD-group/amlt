import pytest
from ase import Atoms
from amlt import contour_exploration
import numpy as np
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms


def test_kappa1():
    
    radius = 2.6 
    atoms = Atoms('AlAl', positions=[[0, 0, 0],[radius,0,0]])
    atoms.set_constraint(FixAtoms(indices=[0]))
    atoms.center(vacuum=10)
    atoms.calc = EMT()

    atoms.set_momenta([[0,0,0],[0,10,0]])

    atoms.set_calculator(EMT())
    initial_energy = atoms.get_potential_energy()

    print("Target Radius (1/kappa) {: .6f} Ang".format( radius))

    name = 'test_kappa1'
    traj_name = name+'.traj'
    log_name = name+'.log'


    dyn = contour_exploration(atoms,
                    maxstep = 1.5,
                    parallel_drift = 0.0,
                    remove_translation  = False,
                    force_parallel_step_scale = None,
                    energy_target = initial_energy,
                    angle_limit = 30,
                    use_tangent_curvature= False,
                    #trajectory = traj_name,
                    #logfile = log_name,
                    ) 

    for i in range(5):
        dyn.run(30)
        
        print('Radius (1/kappa) {: .6f} Ang'.format( 1/dyn.kappa))
        assert 0 == pytest.approx(radius - 1/dyn.kappa, abs=5e-3)
        


def test_kappa2():
    
    pair_distance = 2.6 
    radius = pair_distance*np.sqrt(2)/2
    atoms = Atoms('AlAl', positions=[[-pair_distance/2, 0, 0],[pair_distance/2,0,0]])
    atoms.center(vacuum=10)
    atoms.calc = EMT()

    atoms.set_momenta([[0,-10,0],[0,10,0]])

    atoms.set_calculator(EMT())
    initial_energy = atoms.get_potential_energy()

    print("Target Radius (1/kappa) {: .6f} Ang".format( radius))

    name = 'test_kappa2'
    traj_name = name+'.traj'
    log_name = name+'.log'


    dyn = contour_exploration(atoms,
                    maxstep = 1.0,
                    parallel_drift = 0.0,
                    remove_translation  = False,
                    force_parallel_step_scale = None,
                    energy_target = initial_energy,
                    angle_limit = 15,
                    use_tangent_curvature= False,
                    #trajectory = traj_name,
                    #logfile = log_name,
                    ) 

    for i in range(5):
        dyn.run(30)
        
        print('Radius (1/kappa) {: .6f} Ang'.format( 1/dyn.kappa))
        assert 0 == pytest.approx(radius - 1/dyn.kappa, abs=1e-3)


if __name__ == '__main__':
    test_kappa1()
    test_kappa2()

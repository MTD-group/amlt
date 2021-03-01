import pytest
from ase import Atoms
from amlt import contour_exploration
import numpy as np
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms

from ase import io



def test_momentum_conservation1():
    
    pair_distance = 2.6 
    atoms = Atoms('AlAl', positions=[[0,-pair_distance/2, 0],[0,pair_distance/2,0]])
    atoms.center(vacuum=10)
    atoms.calc = EMT()

    atoms.set_momenta([[1,0,0],[0,0,0]])
    momentum_target = atoms.get_momenta().sum(axis = 0)
    
    atoms.set_calculator(EMT())
    initial_energy = atoms.get_potential_energy()

    print("Target Momentum: [{: .6f} {: .6f} {: .6f}]".format( *tuple(momentum_target)))

    name = 'test_momentum_conservation1'
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
                    trajectory = traj_name,
                    #logfile = log_name,
                    ) 

    for i in range(5):
        dyn.run(30)
        momentum = atoms.get_momenta().sum(axis = 0)
        print("Momentum: [{: .6f} {: .6f} {: .6f}]".format( *tuple(momentum) ))
        assert 0 == pytest.approx(np.linalg.norm(momentum_target - momentum), abs=1e-0)



def test_momentum_conservation2():
    
    d = 1.9
    atoms = Atoms('Al4', positions=[[d,0, 0],[0,d,0],[-d,0,0], [0,-d,0]])
    atoms.center(vacuum=10)
    atoms.calc = EMT()

    #atoms.set_momenta([[1,0,0],[0,0,0],[0,0,0],[0,0,0]])
    atoms.set_momenta([[0,1,0],[-1,0,0],[0,-1,0],[1,0,0]])
    atoms.set_momenta (atoms.get_momenta()*1 + np.array([0,0,4]))
    momentum_target = atoms.get_momenta().sum(axis = 0)
    
    atoms.set_calculator(EMT())
    initial_energy = atoms.get_potential_energy()
    
    #io.write('debug.traj', atoms)
    
    print("Target Momentum: [{: .6f} {: .6f} {: .6f}]".format( *tuple(momentum_target)))

    name = 'test_momentum_conservation2'
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
                    trajectory = traj_name,
                    #logfile = log_name,
                    ) 

    for i in range(10):
        dyn.run(30)
        momentum = atoms.get_momenta().sum(axis = 0)
        print("Momentum: [{: .6f} {: .6f} {: .6f}]".format( *tuple(momentum) ))
        assert 0 == pytest.approx(np.linalg.norm(momentum_target - momentum), abs=1e2)


if __name__ == '__main__':
    test_momentum_conservation1()
    test_momentum_conservation2()

import pytest
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from amlt import contour_exploration
import numpy as np
from ase.calculators.emt import EMT


def test_potentiostat():
    
    size = 2
    atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                              symbol='Al',
                              size=(size, size, size),
                              pbc=True)

    atoms.set_calculator(EMT())
    E0 = atoms.get_potential_energy()

    atoms.rattle(stdev=0.18 , seed = 312)
    
    md_temp = 300
    rng = np.random.RandomState(60622)
    MaxwellBoltzmannDistribution(atoms, temperature_K=md_temp, rng=rng)


    initial_energy = atoms.get_potential_energy()

    print("Energy Above Ground State: {: .4f} eV/atoms".format(
        (initial_energy-E0)/len(atoms)))

    name = 'test_potentiostat'
    traj_name = name+'.traj'
    log_name = name+'.log'


    dyn = contour_exploration(atoms,
                    maxstep = 1.0,
                    parallel_drift = 0.05,
                    remove_translation  = True,
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

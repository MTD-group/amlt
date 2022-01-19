

import numpy as np
rng = np.random.default_rng(12345)

import os
from ase import io

from amlt.polymorpher2 import Polymorpher2

name1 = 'TlInS2_mp-632539_primitive.cif'
name2 = 'Al2.OUTCAR'
name_out1 = 'single_option.traj'
name_out2 = 'two_options.traj'


unit_cell1 = io.read(os.path.abspath('../test_structures/'+name1))
unit_cell2 = io.read(os.path.abspath('../test_structures/'+name2))


#########
pmorpher1 = Polymorpher2(unit_cell1,  rcut = 5.0, rng=rng)
traj = io.Trajectory(name_out1, 'w')
traj.write(unit_cell1)
for i in range(10):
    traj.write(pmorpher1())
traj.close()



#########
pmorpher2 = Polymorpher2([unit_cell1,unit_cell2],  rcut = 6.0, rng=rng)
traj = io.Trajectory(name_out2, 'w')
traj.write(unit_cell1)
traj.write(unit_cell2)
for i in range(4):
    traj.write(pmorpher2())
traj.close()



import numpy as np
import os
from ase import io

from amlt.polymorpher2 import Polymorpher2

name = 'TlInS2_mp-632539_primitive.cif'
name_out = name.replace('.cif', '.traj' )

unit_cell = io.read(os.path.abspath('../test_structures/'+name))

rng = np.random.default_rng(12345)

pmorpher = Polymorpher2(unit_cell,  rcut = 5.0, rng=rng)



traj = io.Trajectory(name_out, 'w')
traj.write(unit_cell)
for i in range(10):
    traj.write(pmorpher())
traj.close()

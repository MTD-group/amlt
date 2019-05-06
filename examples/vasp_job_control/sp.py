import os
os.environ['VASP_PP_PATH'] = '/global/homes/m/mjwaters/vasp_pseudos'

from ase import io
atoms = io.read('POSCAR.initial', format = 'vasp')

from amlt import kgrid_from_cell_volume
kpd = 1000
kpts = kgrid_from_cell_volume(atoms,kpd)
if sum(kpts) == 3:
	#os.environ['VASP_COMMAND'] ="srun -c4 --cpu_bind=cores vasp_gam"
	os.environ['VASP_COMMAND'] ="srun --cpu_bind=cores vasp_gam"
else:
	#os.environ['VASP_COMMAND'] ="srun -c4 --cpu_bind=cores vasp_std"
	os.environ['VASP_COMMAND'] ="srun --cpu_bind=cores vasp_std"	


maxmove = 0.6

### minimal tags to run
from ase.calculators.vasp import Vasp
calc = Vasp( xc='PBE', setups= 'minimal', kpts = kpts)

### good tags to set
calc.set(prec='Accurate', ediff=1E-8*len(atoms), encut = 450, algo = 'Fast', ispin = 1, nelm = 200, lmaxmix = 4)
calc.set(lreal='Auto')


	
### vasp settings
calc.set(nsw = 0, ediffg = -0.01, isym=0 )
## VTST settings
calc.set(ibrion = 3, potim = 0, iopt = 1, maxmove = maxmove )

calc.calculate(atoms)



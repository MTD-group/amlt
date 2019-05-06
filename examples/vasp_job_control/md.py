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


#####
from numpy import loadtxt
temp = loadtxt('temperature.txt')

### minimal tags to run
from ase.calculators.vasp import Vasp
calc = Vasp( xc='PBE', setups= 'minimal', kpts = kpts)

### good tags to set
calc.set(prec='Accurate', ediff= 1E-8*len(atoms), encut = 450, algo = 'Fast', ispin = 1, nelm = 200, lmaxmix = 4)
calc.set(lreal='Auto')


	
################ molecular dynamics  ###########

### vasp settings
calc.set(isym = 0, ibrion = 0) # turn on MD
calc.set(nsw = 10000, potim = 1.0) # 10000 steps, 2.0 fs step

#calc.set(smass=-3, tebeg = 300 )#  NVE, tebeg sets the intial random velocities
calc.set(smass=0, tebeg = temp[0], teend = temp[1]  )# Nose thermostat, NVT

calc.calculate(atoms)



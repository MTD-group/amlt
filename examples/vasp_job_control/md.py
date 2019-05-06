
from ase import io
from ase.calculators.vasp import Vasp
from numpy import loadtxt

temp = loadtxt('temperature.txt')

atoms = io.read('POSCAR.initial', format = 'vasp')


### minimal tags to run
calc = Vasp( xc='PBE', setups= 'minimal', kpts = [1,1,1])

### good tags to set
calc.set(prec='Accurate', ediff= 1E-8*len(atoms), encut = 450, algo = 'Fast', ispin = 1, nelm = 200, lmaxmix = 4)
calc.set(lreal='Auto')


	
################ molecular dynamics  ###########

### vasp settings
calc.set(isym = 0, ibrion = 0) # turn on MD
calc.set(nsw = 10000, potim = 2.0) # 10000 steps, 2.0 fs step

#calc.set(smass=-3, tebeg = 300 )#  NVE, tebeg sets the intial random velocities
calc.set(smass=0, tebeg = temp[0], teend = temp[1]  )# Nose thermostat, NVT

calc.calculate(atoms)



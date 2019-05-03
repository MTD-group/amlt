
from ase import io
from ase.calculators.vasp import Vasp

#max_move_per_atom = 0.05
#maxmove = max_move_per_atom * len(atoms)
maxmove = 0.6

atoms = io.read('POSCAR.initial', format = 'vasp')


### minimal tags to run
calc = Vasp( xc='PBE', setups= 'minimal', kpts = [1,1,1])

### good tags to set
calc.set(prec='Accurate', ediff=1E-7, encut = 450, algo = 'Fast', ispin = 1, nelm = 200, lmaxmix = 4)
calc.set(lreal='Auto')


	
### vasp settings
calc.set(nsw = 0 )

calc.calculate(atoms)



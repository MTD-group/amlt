

#################


from ase import io




submit = False


rcut = 5.0

name_prefix = 'ZrO2'


jobs = [[3, 'random',      'relax'],
        [4, 'polymorphD3', 'md'   ],
        [0, 'known',       'sp'   ]]# when 'known' is seen, well just enumerate the known structures. 




########## random params ######

elements = [ 'O', 'Zr']

element_radii = {	'O' : 1.226399/2,
			'Zr': 2.246900/2}

hard_radii = {		'O' : 1.226399*0.9/2,
                        'Zr': 2.246900*0.9/2}

def generate_random_composition(elements, no_gas = True):
	from numpy.random import rand
	atomic_fractions = rand(len(elements))
	atomic_fractions = atomic_fractions/atomic_fractions.sum()

	if no_gas:
		while atomic_fractions[0] > 0.70:
			atomic_fractions = rand(len(elements))
			atomic_fractions = atomic_fractions/atomic_fractions.sum()

	return atomic_fractions


random_structure_parameters = dict(
		elements = elements, 
		fill_factor_max = 0.30, # 0.65 is about the max for the hard radii
		fill_factor_min = 0.10,
		composition_generator = generate_random_composition,
		cut_off_radius = rcut,
		element_radii = element_radii,
		hard_radii    = hard_radii)


####### polymorphD3 params

from glob import glob

known_structures = []
known_structure_file_names = sorted(glob('known_structures/*.cif'))
for structure_file in known_structure_file_names:
	known_structures.append(io.read(structure_file))

polymorphD3_parameters = dict(
	atom_distortion = 0.2, 
	lattice_distortion = 0.10, 
	deletion_chance = 0.05, 
	rcut = rcut)



##################

job_command = 'sbatch'

job_script_name = 'job.sbatch'

job_script_template = \
'''#!/bin/bash
#SBATCH  --job-name="%s" 
#SBATCH  --output=stdout.log 
#SBATCH  --qos=regular
#SBATCH  --nodes=1 
#SBATCH  --core-spec=4 
#SBATCH  --tasks-per-node=64 
#SBATCH  --time=01:00:00 
#SBATCH  --constraint=knl 
#SBATCH  --account=m3179 
#SBATCH  --mail-type=ALL 
#SBATCH  --mail-user=michael.j.waters@northwestern.edu 

date
pwd

module load python/3.6-anaconda-5.2
#module load vasp/5.4.4-knl
#module load vasp/20171017-knl
module load vasp-tpc/20170629-knl

export OMP_NUM_THREADS=4


export VASP_COMMAND="srun -c4 --cpu_bind=cores vasp_gam"
export VASP_PP_PATH=$HOME/vasp_pseudos

echo "VASP pseudo path" $VASP_PP_PATH

python ../../%s.py
'''


###############
from amlt import *

vasp_job_maker(name_prefix = name_prefix, jobs = jobs,
	job_command = job_command, job_script_name=job_script_name, job_script_template =job_script_template,
	md_temperature_range = (200, 1800),
	submit = submit,
	random_structure_parameters = random_structure_parameters, 
	known_structures=known_structures, polymorphD3_parameters =  polymorphD3_parameters )


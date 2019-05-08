

n_structures = 200 # half will be trained on

use_tensorflow = False

potential_file = 'Zr_O.amp'


cut_off_radius = 6.5





import os
from ase import io


import numpy as np

from amp.utilities import assign_cores




from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.descriptor.zernike import Zernike
from amp.descriptor.bispectrum import Bispectrum
from amp.descriptor.cutoffs import Cosine, Polynomial
import time


if use_tensorflow:
	from amp.model.tflow import NeuralNetwork
else:
	from amp.model.neuralnetwork import NeuralNetwork
from amp.model.kernelridge import KernelRidge

from amp.model import LossFunction


collect_in_zip = True

struct_types = ['known',
		#'polymorphD3',
		'random']
		
dyn_types = ['md','relax','sp']		

from os import path, chdir, getcwd, system
basedir = getcwd()

from ase import io

total = 0

from glob import glob




assign_cores(cores={'localhost': 8})

from amp.descriptor.gaussian import  make_symmetry_functions

##############################################################
if True:
	elements = ['Zr', 'O']
	
	#cutoff=Polynomial(cut_off_radius,gamma = 4)
	cutoff=Cosine(cut_off_radius)
	
	sigmas = np.logspace(np.log10(0.03), np.log10(0.8), num = 9)
	etas = 1.0/(2.0*sigmas**2)
	Gf = make_symmetry_functions(elements=elements, type='G2',
						etas=etas)

	sigmas = np.logspace(np.log10(0.03), np.log10(0.8), num = 7)
	etas = 1.0/(2.0*sigmas**2)
	zetas  = np.logspace(np.log10(1), np.log10(32), num = 4)
	Gf += make_symmetry_functions(elements=elements, type='G4',
						etas=etas,
						zetas=zetas,
						gammas=[1.0, -1.0]) #lambda on the AMP docs

	G = {'Zr': Gf,
		'O': Gf}
	# NN
	#layers = 3
	#nodes = 10
	#layer_config = layers*[nodes]
	layer_config = (16, 8, 4)
	if use_tensorflow:
		my_model = NeuralNetwork(  hiddenlayers={"Zr":layer_config, "O":layer_config}, force_coefficient=None,  )
	else:
		my_model = NeuralNetwork(  hiddenlayers={"Zr":layer_config, "O":layer_config},  checkpoints = 20)


if os.path.isfile(potential_file ):
	print("Loading", potential_file )
	MLIP = Amp.load(potential_file )
else:
	print("Training New Potential")
	MLIP = Amp(descriptor=Gaussian(Gs=G,
				cutoff=cutoff),
				model=my_model)



from amlt import super_cell_if_needed

image_list = []

time1 = time.time()
for struct_type in struct_types:
	for dyn_type in dyn_types:	
		top_direct = ('%s_%s/')%(struct_type, dyn_type)
		if path.isdir(top_direct):
			print(top_direct)
			sub_total = 0
			sub_dirs = sorted(glob(top_direct+'*/'))
			for sub_dir in sub_dirs:
				name = sub_dir.split('/')[-2]
				if name.isdigit():
					
					traj = io.Trajectory(filename = sub_dir  + 'images.traj', mode='r')

					print (sub_dir.ljust(22)+ ('atoms %i'%len(traj[0])).ljust(12) +'images', len(traj))
					
					for image in traj:
						im = super_cell_if_needed(image, cut_off_radius, verbose=False)						
						image_list.append(im)
						
						
					
					sub_total += len(traj)
			print('sub_total: %i \n'% sub_total)

			total+=	sub_total
print('Total number of images:', total)
time2 = time.time()
print('Time for file parsing is: {}'.format(time2-time1))

##################################################


# quick train with energy?
MLIP.model.lossfunction = LossFunction(convergence={'energy_rmse': 0.001}, force_coefficient=None)
MLIP.train(images=image_list[0:200])

if False:
	# slow train with energy and forces?
	MLIP.model.lossfunction = LossFunction( energy_coefficient=1.0, force_coefficient=0.2,
								convergence={'energy_rmse': 0.1, 'force_rmse': 0.6})
	MLIP.train(images=image_list)

MLIP.save(potential_file, overwrite=True)



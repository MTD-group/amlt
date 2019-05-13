



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




struct_types = [#'known',
		#'polymorphD3']
		'random']
		
dyn_types = ['md','relax','sp']		

from os import path, chdir, getcwd, system
basedir = getcwd()

from ase import io
from glob import glob






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
	#layers = 1
	#nodes = 10
	#layer_config = layers*[nodes]
	layer_config = (16, 8, 4)
	if use_tensorflow:
		my_model = NeuralNetwork(  hiddenlayers={"Zr":layer_config, "O":layer_config}, force_coefficient=None,  )
	else:
		my_model = NeuralNetwork(  hiddenlayers={"Zr":layer_config, "O":layer_config},  checkpoints = 5)


if os.path.isfile(potential_file ):
	print("Loading", potential_file )
	MLIP = Amp.load(potential_file )
else:
	print("Training New Potential")
	MLIP = Amp(descriptor=Gaussian(Gs=G,
				cutoff=cutoff),
				model=my_model)



############### parsing images ##########

from amlt import super_cell_if_needed

image_list = []
total=0
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
					
					if int(name)%2 == 0:
						traj = io.Trajectory(filename = sub_dir  + 'images_supercell.traj', mode='r')
						#traj_super = io.Trajectory(filename = sub_dir  + 'images_supercell.traj', mode='w')
						
						subsub_total = 0
					
						#for image in traj:
						for image_index in range(len(traj)):
							if image_index >= 5 and image_index%1 == 0:
								image = traj[image_index]
							
								verbose = False
								if image_index == 0: verbose = True
							
								im = super_cell_if_needed(image, cut_off_radius, verbose=verbose)						
								#traj_super.write(im)
								image_list.append(im)
						
								subsub_total += 1
								
						print (sub_dir.ljust(22)+ ('atoms %i'%len(traj[0])).ljust(12) +'images loaded %i/%i'%(subsub_total, len(traj)))
						sub_total += subsub_total
						#traj_super.close()
						traj.close()
			print('sub_total: %i \n'% sub_total)

			total+=	sub_total
print('Total Number of Training Images:', total)
time2 = time.time()
print('Time for file parsing is: {}'.format(time2-time1))

##################################################

from random import shuffle
# quick train with energy?
MLIP.model.lossfunction = LossFunction(convergence={'energy_rmse': 0.001}, force_coefficient=None)
shuffle(image_list)
#assign_cores(cores={'localhost': 16})
MLIP.train(images=image_list)

if False:
	# slow train with energy and forces?
	MLIP.model.lossfunction = LossFunction( energy_coefficient=1.0, force_coefficient=0.2,
								convergence={'energy_rmse': 0.1, 'force_rmse': 0.6})
	MLIP.train(images=image_list)

MLIP.save(potential_file, overwrite=True)



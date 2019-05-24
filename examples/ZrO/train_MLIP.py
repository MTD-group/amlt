
potential_file = 'Zr_O.amp'

elements = ['Zr', 'O']
cut_off_radius = 6.5

overfit_penalty = 0.2

checkpoint_interval = 4
use_tensorflow = False
cores = 8


log_file = 'amp-train.log'

########### loading AMP file if present ########################################
from amp import Amp
import os


if os.path.isfile(potential_file ):
    print("Loading", potential_file )
    MLIP = Amp.load(potential_file )

else: ############ Setting Up New Descriptors and MLIP ################

    from amp.descriptor.gaussian import Gaussian
    from amp.descriptor.zernike import Zernike
    from amp.descriptor.bispectrum import Bispectrum
    from amp.descriptor.cutoffs import Cosine, Polynomial
    from amp.descriptor.gaussian import  make_symmetry_functions
    from amp.model.kernelridge import KernelRidge
    import numpy as np
    
    print("Training New Potential")

    #cutoff=Polynomial(cut_off_radius,gamma = 4)
    cutoff=Cosine(cut_off_radius)
    
    sigmas = np.logspace(np.log10(0.03), np.log10(0.8), num = 9)
    etas = 1.0/(2.0*sigmas**2)
    Gf = make_symmetry_functions(elements=elements, type='G2',
                        etas=etas)

    sigmas = np.logspace(np.log10(0.03), np.log10(0.8), num = 7)
    zetas  = np.logspace(np.log10(1.0),  np.log10(32),  num = 4)
    etas = 1.0/(2.0*sigmas**2)
    Gf += make_symmetry_functions(elements=elements, type='G4',
                                etas=etas,
                                zetas=zetas,
                                gammas=[1.0, -1.0]) #lambda on the AMP docs

    layer_config = (16, 8, 4)

    G = {}
    hiddenlayers = {}
    for el in elements:
        G[el] =  Gf
        hiddenlayers[el] = layer_config

    if use_tensorflow:
        from amp.model.tflow import NeuralNetwork
        my_model = NeuralNetwork(  hiddenlayers=hiddenlayers, force_coefficient=None  )
    else:
        from amp.model.neuralnetwork import NeuralNetwork
        my_model = NeuralNetwork(  hiddenlayers=hiddenlayers)

    MLIP = Amp(descriptor=Gaussian(Gs=G,
                cutoff=cutoff),
                model=my_model)


MLIP.model.checkpoints = checkpoint_interval
MLIP.cores['localhost'] = cores
from amp.utilities import Logger, make_filename
MLIP._log = Logger(make_filename('',log_file))

############### parsing images ##########
from ase import io
from glob import glob
import time

struct_types = [ 'known', 'polymorphD3', 'random' ]
dyn_types    = [ 'md','relax','sp']        

image_list = []
total = 0
time1 = time.time()
for struct_type in struct_types:
    for dyn_type in dyn_types:    
        top_direct = ('%s_%s/')%(struct_type, dyn_type)
        if os.path.isdir(top_direct):
                print(top_direct)
                sub_total = 0
                sub_dirs = sorted(glob(top_direct+'*/'))
                for sub_dir in sub_dirs:
                    name = sub_dir.split('/')[-2]
                    if name.isdigit():
                        traj = io.Trajectory(filename = sub_dir  + 'images_supercell.traj', mode='r')
                        subsub_total = 0
                        for image_index in range(len(traj)):
                            if int(name)%2 == 0:
                            #if image_index >= 0 and image_index%1 == 0:
                                image = traj[image_index]
                                image_list.append(image)
                                subsub_total += 1
                            elif struct_type == 'known':
                                image = traj[image_index]
                                image_list.append(image)
                                subsub_total += 1
                            
                        print (sub_dir.ljust(22)+ ('atoms %i'%len(traj[0])).ljust(12) +'images loaded %i/%i'%(subsub_total, len(traj)))
                        sub_total += subsub_total
                        traj.close()
                # subtotal for this struct+dyn_type
                print('sub_total: %i \n'% sub_total)
                total+=    sub_total


print('Total Number of Training Images:', total)
print('Time for file parsing is: {}'.format(time.time() - time1))

##################################################
from amp.model import LossFunction
from random import shuffle
def train_model(energy_coeff=1.0, force_training=False, force_coeff=0.2):
    shuffle(image_list) # "Load balancing"
    if force_training:
        # slow train with energy and forces?
        MLIP.model.lossfunction = LossFunction(
                energy_coefficient=energy_coeff, 
                force_coefficient=force_coeff,
                convergence={'energy_rmse': 0.1, 'force_rmse': 0.6})
    else:
        MLIP.model.lossfunction = LossFunction(
             convergence={'energy_rmse': 0.001, 'overfit': overfit_penalty},
             force_coefficient=None, overfit = overfit_penalty)
    
    MLIP.train(images=image_list) 
    MLIP.save(potential_file, overwrite=True)
train_model()


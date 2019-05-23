


from ase import io
#from ase.calculators.emt import EMT
#from ase.build import fcc110
#from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
#from ase.md import VelocityVerlet

from amp import Amp
from amlt import super_cell_if_needed


import time

cut_off_radius = rcut = 6.5


#MLIP = Amp.load('Zr_O.amp')
MLIP = Amp.load('amp-checkpoint.amp')

struct_types = [
        'known',
        'polymorphD3',
        'random']

dyn_types = ['md','relax','sp']        

from os import path, getcwd
basedir = getcwd()
from glob import glob


############
#from amlt import super_cell_if_needed

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
                    

                    traj = io.Trajectory(filename = sub_dir  + 'images.traj', mode='r')


    
                    image_index = -1
                    image = traj[image_index]

                    verbose = False
                    

#                    im = super_cell_if_needed(image, cut_off_radius, verbose=verbose)                        
                    im = image
                

                    print (sub_dir.ljust(22)+ ('atoms %i'%len(traj[0])).ljust(12) +'images loaded %i/%i'%(1, len(traj)))
                    DFT_e = image.get_potential_energy() / len(image)
                    NN_e  = MLIP.get_potential_energy(im) / len(im)
                    
                    print("DFT %f NN %f"%(DFT_e,NN_e))

                    traj.close()
            



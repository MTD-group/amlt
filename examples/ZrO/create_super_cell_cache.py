

cut_off_radius = 6.5


############### parsing images ##########
from amlt import super_cell_if_needed
from ase import io
from glob import glob
import os
import time
struct_types = ['known',
				'polymorphD3',
				'random']
		
dyn_types = ['md','relax','sp']		


#image_list = []
total=0
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
					
					traj = io.Trajectory(filename = sub_dir  + 'images.traj', mode='r')
					traj_super = io.Trajectory(filename = sub_dir  + 'images_supercell.traj', mode='w')
					
					subsub_total = 0
				
					for image_index in range(len(traj)):
						image = traj[image_index]
					
						verbose = False
						if image_index == 0: verbose = True
					
						im = super_cell_if_needed(image, cut_off_radius, verbose=verbose)						
						traj_super.write(im)
						
				
						subsub_total += 1
							
					print (sub_dir.ljust(22)+ ('atoms %i'%len(traj[0])).ljust(12) +'images loaded %i/%i'%(subsub_total, len(traj)))
					sub_total += subsub_total
					traj_super.close()
					traj.close()
			print('sub_total: %i \n'% sub_total)

			total+=	sub_total
print('Total Number of Training Images:', total)
time2 = time.time()
print('Time for file parsing is: {}'.format(time2-time1))






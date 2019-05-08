
collect_in_zip = True

struct_types = ['known',
		'polymorphD3',
		'random']
		
dyn_types = ['md','relax','sp']		

from os import path, chdir, getcwd, system
basedir = getcwd()

from ase import io

total = 0

from glob import glob

for struct_type in struct_types:
	for dyn_type in dyn_types:
		
		#chdir(basedir)
		
		top_direct = ('%s_%s/')%(struct_type, dyn_type)
		if path.isdir(top_direct):
			#chdir(top_direct)
			print(top_direct)
			sub_total = 0
			sub_dirs = sorted(glob(top_direct+'*/'))
			for sub_dir in sub_dirs:
				name = sub_dir.split('/')[-2]
				#name =
				#print(name)
				if name.isdigit():
					
					traj = io.Trajectory(filename = sub_dir  + 'images.traj', mode='r')

					print (sub_dir.ljust(22)+ ('atoms %i'%len(traj[0])).ljust(12) +'images', len(traj))
					if collect_in_zip:
						system('zip -r  Zr-O_amlt_data_set.zip ' + sub_dir  + 'images.traj')
					
					sub_total += len(traj)
			print('sub_total: %i \n'% sub_total)

			total+=	sub_total
print('Total:', total)


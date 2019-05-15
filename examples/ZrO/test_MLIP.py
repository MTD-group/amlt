



from ase import io
import os

from amp import Amp
from numpy import array, mean, sqrt, floor, arange, linspace, savetxt, percentile

from amlt import super_cell_if_needed


n_structures = 50

rcut = 6.5
use_uc = False

MLIP = Amp.load('Zr_O.amp')


actual_energy_train =[]
MLIP_energy_train =  []
actual_energy_test =[]
MLIP_energy_test =  []

train_images = []
test_images  = []

for structure_number in range(n_structures):
    struct_dir = os.path.join('random_md', str(structure_number))
    print (struct_dir)
    filepath = os.path.join(struct_dir, 'images_supercell.traj')
    traj =  io.Trajectory(filename=filepath, mode='r')	


    for atoms in traj:
        image = super_cell_if_needed(atoms, rcut =rcut)
        if structure_number%2==0:		
            actual_energy_train.append( image.get_potential_energy()/len(image))
            MLIP_energy_train.append( MLIP.get_potential_energy(image)/len(image))
            train_images.append(image.copy())
        else:
            actual_energy_test.append( image.get_potential_energy()/len(image))
            MLIP_energy_test.append( MLIP.get_potential_energy(image)/len(image))
            test_images.append(image.copy())

    traj.close()



############

savetxt('train_energies.txt', array([actual_energy_train, MLIP_energy_train]).T)
savetxt('test_energies.txt', array([actual_energy_test, MLIP_energy_test]).T)

###########
if use_uc:
	images = io.read('unit_cell/'+'OUTCAR', index = ':')
	for atoms in images:
		image = super_cell_if_needed(atoms, rcut =rcut)

		actual = image.get_potential_energy()/len(image)
		MLIP   =  MLIP.get_potential_energy(image)/len(image)


		actual_energy_test.append(actual)
		MLIP_energy_test.append(  MLIP)
		test_images.append(image)

#chdir(basedir)


from amp.model.neuralnetwork import NodePlot
nodeplot = NodePlot(MLIP)
nodeplot.plot(train_images, filename='node_plot_train.pdf')
#nodeplot.plot(test_images, filename='node_plot_test.pdf') # why does this crash?

############# now we make a comparison plot #############

from matplotlib import pyplot

fig, ax = pyplot.subplots(figsize = (3.2,2.5), dpi = 200)

ax.plot(actual_energy_test, MLIP_energy_test, 'r.', label = 'Test', markersize = 3)
ax.plot(actual_energy_train, MLIP_energy_train, 'b.', label = 'Train', markersize = 3)
if use_uc: ax.text(c_ZrO2_actual, c_ZrO2_MLIP , 'c-ZrO$_{2}$', fontsize = 8)


max_e  = max( max(ax.get_xlim()) ,  max(ax.get_ylim()) )
min_e  = min( min(ax.get_xlim()) ,  min(ax.get_ylim()) )

ax.plot([min_e, max_e], [min_e, max_e], color = 'k', linestyle = '--', zorder = 0)

ax.set_xlabel('Actual Energy (eV/atom)')
ax.set_ylabel('MLIP Energy (eV/atom)')

ax.legend()
ax.minorticks_on()

fig.tight_layout(pad = 0.1)
fig.savefig('parityplot.png')
fig.savefig('parityplot.pdf')










fig, axes = pyplot.subplots(nrows = 2, ncols = 2, figsize = (3.2*2,2.5*2), dpi = 200)


train_error = array(MLIP_energy_train)-array(actual_energy_train)
test_error =  array(MLIP_energy_test)-array(actual_energy_test)

#print(train_error, test_error)
axes[0,0].plot(actual_energy_test, test_error, 'r.', label = 'Test', markersize = 3)
axes[0,0].plot(actual_energy_train, train_error, 'b.', label = 'Train', markersize = 3)
axes[0,0].axhline(0, color = 'grey', linestyle = '--')
axes[0,0].set_xlabel('Actual Energy (eV/atom)')
axes[0,0].set_ylabel('MLIP Error (eV/atom)')
axes[0,0].legend()
axes[0,0].minorticks_on()

if use_uc: axes[0,0].text(c_ZrO2_actual, c_ZrO2_MLIP - c_ZrO2_actual, 'c-ZrO$_{2}$', fontsize = 8)


################
def nice_bins(errors, bin_step_scale = 0.2):


	rms = sqrt(mean(errors**2))
	bin_step = bin_step_scale*rms
	bin_min = bin_step * floor(errors.min()/bin_step )
	bin_max = errors.max()+bin_step
	bins = arange(bin_min,bin_max , bin_step  )
	if len(bins)< 10:
		bins = linspace(errors.min(), errors.max(), 10)
	print(bins[-1]/bin_step, errors.max()/bin_step )
	return bins



def nice_bins_percentile(errors, bins_in_window = 30, window_percent = 90):

	window = percentile(errors, [50 - window_percent/2.0, 50 + window_percent/2.0])
	bin_step = (window.max()-window.min())/bins_in_window
	#rms = sqrt(mean(errors**2))
	#bin_step = bin_step_scale*rms
	bin_min = bin_step * floor(errors.min()/bin_step )
	bin_max = errors.max()+bin_step

	bins = arange(bin_min,bin_max , bin_step  )
	if len(bins)< 10:
		bins = linspace(errors.min(), errors.max(), 10)
	print(bins[-1]/bin_step, errors.max()/bin_step )
	return bins

def auto_limits(errors, scale = 3, window_percent = 90):
	p = percentile(errors, [50 - window_percent/2.0, 50 , 50 + window_percent/2.0])
	return ( p[1]+scale*(p[0]-p[1]) , p[1]+scale*(p[2]-p[1]))



#############################3
test_rms = sqrt(mean(test_error**2))
train_rms = sqrt(mean(train_error**2))

axes[0,1].hist(test_error,  bins = nice_bins_percentile(test_error, bin_step_scale = 0.05),  
				orientation="horizontal", color = 'r', label = 'RMS: %f'%test_rms)
axes[0,1].hist(train_error, bins = nice_bins_percentile(train_error), 
				orientation="horizontal", color = 'b', label = 'RMS: %f'%train_rms)

axes[0,1].set_ylim(auto_limits(test_error) )
axes[0,1].legend()

##############


axes[1,0].plot(actual_energy_test, test_error, 'r.', label = 'Test', markersize = 3)
axes[1,0].plot(actual_energy_train, train_error, 'b.', label = 'Train', markersize = 3)
axes[1,0].axhline(0, color = 'grey', linestyle = '--')
axes[1,0].set_xlabel('Actual Energy (eV/atom)')
axes[1,0].set_ylabel('MLIP Error (eV/atom)')
axes[1,0].legend()
axes[1,0].minorticks_on()

data_range = train_error.max() - train_error.min()
delta = 0.05 * data_range
axes[1,0].set_ylim(train_error.min()-delta , train_error.max()+delta)


##############

axes[1,1].hist(train_error,bins  = nice_bins_percentile(train_error), orientation="horizontal", color = 'b', label = 'RMS: %f'%train_rms)
axes[1,1].set_ylim(auto_limits(train_error))


fig.tight_layout(pad = 0.1)
fig.savefig('error_plot.png')
fig.savefig('error_plot.pdf')







pyplot.show()
###############################################################################

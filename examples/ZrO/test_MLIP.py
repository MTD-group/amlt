

use_cache = False
make_node_plots  = False
MLIP_file_name = 'Zr_O.amp'
rcut = 6.5

log_file = 'amp-test.log'
#############################################################
from ase import io
from glob import glob
import os

struct_types = [
        'known',
        'polymorphD3',
        'random']



dyn_types = ['md','relax','sp']



actual_energy_train =[]
MLIP_energy_train =  []
actual_energy_test =[]
MLIP_energy_test =  []

train_images = []
test_images  = []


total=0
#time1 = time.time()
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
                                    #print()
                                    traj = io.Trajectory(filename = sub_dir  + 'images_supercell.traj', mode='r')
                                    subsub_total = 0
                                    #for image in traj:
                                    for image_index in range(len(traj)):
                                        image = traj[image_index]
                                        if int(name)%2 == 0:
                                        #if image_index >= 0 and image_index%1 == 0:
                                           train_images.append(image)
                                        elif struct_type == 'known':
                                           train_images.append(image)
                                        else:
                                           test_images.append(image)
                                        subsub_total += 1

                                    print (sub_dir.ljust(22)+ ('atoms %i'%len(traj[0])).ljust(12) +'images loaded %i/%i'%(subsub_total, len(traj)))
                                    sub_total += subsub_total
                                    #traj_super.close()
                                    traj.close()
                        print('sub_total: %i \n'% sub_total)

                        total+= sub_total
print('Total Number of Images:', total)

#############################################################

if make_node_plots:
    from amp.model.neuralnetwork import NodePlot
    nodeplot1 = NodePlot(MLIP)
    print('NodePlot Test Images')
    nodeplot1.plot(test_images, filename='node_plot_test.pdf') 
    
    nodeplot2 = NodePlot(MLIP)
    print('NodePlot Training Images')
    nodeplot2.plot(train_images, filename='node_plot_train.pdf')


#############################################################
import numpy as np
if use_cache: 
    actual_energy_train, MLIP_energy_train = np.loadtxt('train_energies.txt', usecols = [0,1], unpack =True)
    actual_energy_test,  MLIP_energy_test  = np.loadtxt('test_energies.txt',  usecols = [0,1], unpack =True)
else:
    from amp import Amp
    from amp.utilities import Logger, make_filename
    MLIP = Amp.load(MLIP_file_name)
    MLIP._log = Logger(make_filename('',log_file))
    
    def MLIP_eval(image):
        return MLIP.get_potential_energy(image)/len(image)

    from multiprocessing import Pool, cpu_count
    my_pool = Pool(cpu_count())
    
    print('Training Images...')
    MLIP_energy_train = my_pool.map(MLIP_eval, train_images)
    for i in range(len(train_images)):
        image = train_images[i]
        #if i%100==0: print(i,'/', len(train_images))
        actual_energy_train.append( image.get_potential_energy()/len(image))
        #MLIP_energy_train.append( MLIP.get_potential_energy(image)/len(image))
    
    print('Test Images...')
    MLIP_energy_test = my_pool.map(MLIP_eval, test_images)
    for i in range(len(test_images)):
        image = test_images[i]
        #if i%100==0: print(i,'/', len(test_images))
        actual_energy_test.append( image.get_potential_energy()/len(image))
        #MLIP_energy_test.append( MLIP.get_potential_energy(image)/len(image))
   

    np.savetxt('train_energies.txt', np.array([actual_energy_train, MLIP_energy_train]).T)
    np.savetxt('test_energies.txt',  np.array([actual_energy_test,  MLIP_energy_test ]).T)



############# now we make a comparison plot #############

from matplotlib import pyplot

fig, ax = pyplot.subplots(figsize = (3.2,2.5), dpi = 200)

ax.plot(actual_energy_test, MLIP_energy_test, 'r.', label = 'Test', markersize = 3)
ax.plot(actual_energy_train, MLIP_energy_train, 'b.', label = 'Train', markersize = 3)
#if use_uc: ax.text(c_ZrO2_actual, c_ZrO2_MLIP , 'c-ZrO$_{2}$', fontsize = 8)


max_e  = max( max(ax.get_xlim()) ,  max(ax.get_ylim()) )
min_e  = min( min(ax.get_xlim()) ,  min(ax.get_ylim()) )

ax.plot([min_e, max_e], [min_e, max_e], color = 'k', linestyle = '--', zorder = 0)

ax.set_xlabel('DFT Energy (eV/atom)')
ax.set_ylabel('MLIP Energy (eV/atom)')

ax.legend()
ax.minorticks_on()

fig.tight_layout(pad = 0.1)
fig.savefig('parityplot.png', transparent = True, dpi = 600)
fig.savefig('parityplot.pdf', transparent = True, dpi = 600)










fig, axes = pyplot.subplots(nrows = 2, ncols = 2, figsize = (3.2*2,2.5*2), dpi = 200)


train_error = np.array(MLIP_energy_train) - np.array(actual_energy_train)
test_error =  np.array(MLIP_energy_test)  - np.array(actual_energy_test)

#print(train_error, test_error)
axes[0,0].plot(actual_energy_test, test_error, 'r.', label = 'Test', markersize = 3)
axes[0,0].plot(actual_energy_train, train_error, 'b.', label = 'Train', markersize = 3)
axes[0,0].axhline(0, color = 'grey', linestyle = '--')
axes[0,0].set_xlabel('DFT Energy (eV/atom)')
axes[0,0].set_ylabel('MLIP Error (eV/atom)')
axes[0,0].legend()
axes[0,0].minorticks_on()




################
def nice_bins(errors, bin_step_scale = 0.2):


    rms = np.sqrt(np.mean(errors**2))
    bin_step = bin_step_scale*rms
    bin_min = bin_step * np.floor(errors.min()/bin_step )
    bin_max = errors.max()+bin_step
    bins = np.arange(bin_min,bin_max , bin_step  )
    if len(bins)< 10:
        bins = np.linspace(errors.min(), errors.max(), 10)
    #print(bins[-1]/bin_step, errors.max()/bin_step )
    return bins



def nice_bins_percentile(errors, bins_in_window = 30, window_percent = 90):

    window = np.percentile(errors, [50 - window_percent/2.0, 50 + window_percent/2.0])
    bin_step = (window.max()-window.min())/bins_in_window
    #rms = sqrt(mean(errors**2))
    #bin_step = bin_step_scale*rms
    bin_min = bin_step * np.floor(errors.min()/bin_step )
    bin_max = errors.max()+bin_step

    bins = np.arange(bin_min,bin_max , bin_step  )
    if len(bins)< 10:
        bins = np.linspace(errors.min(), errors.max(), 10)
    print(bins[-1]/bin_step, errors.max()/bin_step )
    return bins

def auto_limits(errors, scale = 3, window_percent = 90):
    p = np.percentile(errors, [50 - window_percent/2.0, 50 , 50 + window_percent/2.0])
    return ( p[1]+scale*(p[0]-p[1]) , p[1]+scale*(p[2]-p[1]))



#############################3
test_rms =  np.sqrt(np.mean( test_error**2))
train_rms = np.sqrt(np.mean(train_error**2))

axes[0,1].hist(test_error,  bins = nice_bins_percentile(test_error),   density = True,
                orientation="horizontal", color = 'r', label = 'RMS: %f'%test_rms)
axes[0,1].hist(train_error, bins = nice_bins_percentile(train_error), density = True,
                orientation="horizontal", color = 'b', label = 'RMS: %f'%train_rms)

axes[0,1].set_ylim(auto_limits(test_error) )
axes[0,1].legend()

##############


axes[1,0].plot(actual_energy_test, test_error, 'r.', label = 'Test', markersize = 3)
axes[1,0].plot(actual_energy_train, train_error, 'b.', label = 'Train', markersize = 3)
axes[1,0].axhline(0, color = 'grey', linestyle = '--')
axes[1,0].set_xlabel('DFT Energy (eV/atom)')
axes[1,0].set_ylabel('MLIP Error (eV/atom)')
axes[1,0].legend()
axes[1,0].minorticks_on()

data_range = train_error.max() - train_error.min()
delta = 0.05 * data_range
axes[1,0].set_ylim(train_error.min()-delta , train_error.max()+delta)


##############

axes[1,1].hist(train_error,bins  = nice_bins_percentile(train_error), orientation="horizontal", color = 'b',density = True,  label = 'RMS: %f'%train_rms)
axes[1,1].set_ylim(auto_limits(train_error))



fig.tight_layout(pad = 0.2)
fig.subplots_adjust(left= 0.12)
fig.savefig('error_plot.png', transparent = True, dpi = 600)
fig.savefig('error_plot.pdf', transparent = True, dpi = 600)







pyplot.show()
###############################################################################

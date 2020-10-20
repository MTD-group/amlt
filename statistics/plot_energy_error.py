
import numpy as np
from . import struct_colors, dyn_markers, dyn_types, struct_types
from . import read_evaluation_data, get_energy_lists
from . import nice_bins_percentile




def plot_energy_error(ax1, ax2, data_sets, 
                struct_types = struct_types,
                dyn_types = dyn_types,
                struct_colors = struct_colors,
                dyn_markers = dyn_markers,
                bad_data_traj_list = [],
                zorders = None):
    
    if zorders is None:
        zorders = [z for z in range(1,len(data_sets)+1)]
    
    from matplotlib import colors as mcolors
    def get_color(name):
        return np.array(mcolors.to_rgb(name))
            
    
    for struct_type in struct_types:
    
        color = np.array( get_color(struct_colors[struct_type])  )
    
        for di, data_set in enumerate(data_sets):
            #     fname       name    color  lightness  zorder  
            #       0             1      2      3        4
            fname =  data_set[0]
            data_name = data_set[1]
            lightness = data_set[3]
            zorder = zorders[di]
            
            for dyn_type in dyn_types:
                marker = dyn_markers[dyn_type]
                label = data_name+'\n'+ struct_type +':' + dyn_type
                
                image_pairs = read_evaluation_data(filename = fname, struct_types = [struct_type], dyn_types = [dyn_type])
                if len(image_pairs)>0:
                    cache_energy, data_energy = get_energy_lists(image_pairs)
                    ax1.plot(data_energy, cache_energy-data_energy,
                        color = color*lightness,
                        marker = marker,  markersize = 3, linestyle = '',
                        label = label, zorder = zorder)



    ax1.axhline(0, color = 'grey', linestyle = '--', zorder = -1)
    ax1.set_xlabel('DFT Energy (eV/atom)')
    ax1.set_ylabel('MLIP Error (eV/atom)')
    ax1.legend(fontsize = 8, handletextpad = 0.3)
    
    #ylims = auto_limits(test_error)
    #ylims = (-0.3, 0.7)
    #axes[0].set_ylim( ylims)

    ax1.minorticks_on()
    ylims = ax1.get_ylim()

    
    for di, data_set in enumerate(data_sets):

        fname =  data_set[0]
        data_name = data_set[1]
        color = np.array(get_color(data_set[2]))
        lightness = data_set[3]
        zorder = zorders[di]
        
        
        #for dyn_type in dyn_types:
        #    marker = dyn_markers[dyn_type]
        #label = data_name#+'-'+ struct_type +'-' + dyn_type
        
        image_pairs = read_evaluation_data(filename = fname, struct_types = struct_types, dyn_types = dyn_types)
        cache_energy, data_energy = get_energy_lists(image_pairs)
        error = cache_energy - data_energy
        error_rms = np.sqrt(np.mean( error**2))
        
        my_bins =  nice_bins_percentile(error, bins_in_window = 30, window_percent = 90)
        ax2.hist(error, my_bins, density = False,
            orientation="horizontal", color = color*lightness, label = '%s\nRMS: %.1f meV'%(data_name ,error_rms*1000), zorder = zorder )



    ax2.set_ylim(ylims )
    ax2.legend(fontsize = 8, handletextpad = 0.3, borderpad = 0.1)
    ax2.minorticks_on()
    ax2.set_xlabel('# of Samples')
    #pyplot.setp(ax2.get_yticklabels(), visible=False)
    #ax2.yaxis.set_ticklabels([])


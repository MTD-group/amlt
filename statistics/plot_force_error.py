import numpy as np
from . import struct_colors, dyn_markers, dyn_types, struct_types
from . import read_evaluation_data, get_force_list, compute_force_error_list
from . import compute_rms_force_error_by_atom, compute_rms_force_error_by_image
from . import compute_force_norms_by_image, collapse_sub_lists
from . import nice_bins_percentile

def plot_force_error(ax1, ax2, data_sets, 
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
                    #cache_energy, data_energy = get_energy_lists(image_pairs)
                    cache_forces, data_forces = get_force_list(image_pairs)
                    force_error_list =  compute_force_error_list(cache_forces, data_forces)
                    rms_force_error_by_atom = compute_rms_force_error_by_atom(cache_forces, data_forces)
                    rms_force_error_by_image = compute_rms_force_error_by_image(cache_forces, data_forces)/np.sqrt(3)
                    x, y = compute_force_norms_by_image(data_forces)/np.sqrt(3), rms_force_error_by_image
                    ax1.plot(x,y,
                        color = color*lightness,
                        marker = marker,  markersize = 3, linestyle = '',
                        label = label, zorder = zorder)


    ax1.set_xlabel('DFT Force (eV/Å)')
    ax1.set_ylabel('MLIP Force Error (eV/Å)')
    
    ymax = ax1.get_ylim()[1]
    ax1.set_ylim(0,ymax)
    
    xmax = ax1.get_xlim()[1]
    ax1.set_xlim(0,xmax)

    #ax1.axhline(0, color = 'grey', linestyle = '--', zorder = -1)
    ax1.legend(fontsize = 8 , loc = 'upper left',  borderpad = 0.1)
    ax1.minorticks_on()
    ylims = ax1.get_ylim()
    
    
    
    
    
    ###### histogram part
    
    #structs_of_sets = ['random', 'random', 'polymorphD3']
    
    for di, data_set in enumerate(data_sets):


        fname =  data_set[0]
        data_name = data_set[1]
        color = np.array(get_color(data_set[2]))
        lightness = data_set[3]
        zorder = zorders[di]

        #for dyn_type in dyn_types:
        #    marker = dyn_markers[dyn_type]
        #label = data_name#+'-'+ struct_type +'-' + dyn_type
        #label = data_name+'\n'+ struct_type +':' + dyn_type
        
        image_pairs = read_evaluation_data(filename = fname, struct_types = struct_types, dyn_types = dyn_types)
        cache_forces, data_forces = get_force_list(image_pairs)
        force_error_list =  compute_force_error_list(cache_forces, data_forces)
        rms_force_error_by_atom = compute_rms_force_error_by_atom(cache_forces, data_forces)
        rms_force_error_by_image = compute_rms_force_error_by_image(cache_forces, data_forces)/np.sqrt(3)
        
        net_rms_force_error_by_atom =  np.sqrt(np.mean( collapse_sub_lists(rms_force_error_by_atom)**2))
        net_rms_force_error_by_image = np.sqrt(np.mean( (rms_force_error_by_image)**2))
        
        my_bins =  nice_bins_percentile(rms_force_error_by_image, bins_in_window = 30, window_percent = 90)
        ax2.hist(rms_force_error_by_image, my_bins, density = False,
            orientation="horizontal", color = color*lightness, label = '%s\nRMS: %.1f meV/Å'%(data_name ,net_rms_force_error_by_image*1000), zorder = zorder )



    ax2.set_ylim(ylims )
    ax2.legend(fontsize = 8, handletextpad = 0.3, borderpad = 0.1)
    ax2.minorticks_on()
    ax2.set_xlabel('# of Samples')
    #pyplot.setp(ax2.get_yticklabels(), visible=False)
    #ax2.yaxis.set_ticklabels([])

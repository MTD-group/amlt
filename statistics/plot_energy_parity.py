
import numpy as np
from . import struct_colors, dyn_markers, dyn_types, struct_types
from . import read_evaluation_data, get_energy_lists


def plot_energy_parity(ax, data_sets, 
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
            fname =  data_set[0]
            data_name = data_set[1]
            lightness = data_set[3]
            zorder = zorders[di]
            
            for dyn_type in dyn_types:
                marker = dyn_markers[dyn_type]
                label = data_name+'-'+ struct_type +'-' + dyn_type
                
                image_pairs = read_evaluation_data(filename = fname, 
                        struct_types = [struct_type], 
                        dyn_types = [dyn_type],
                        use_forces = False,
                        bad_data_traj_list = bad_data_traj_list)
                if len(image_pairs)>0:
                    cache_energy, data_energy = get_energy_lists(image_pairs)
                    ax.plot(data_energy, cache_energy,
                        color = color*lightness,
                        marker = marker,  markersize = 3, linestyle = '',
                        label = label, zorder = zorder)
    

    ######

    max_e  = max( max(ax.get_xlim()) ,  max(ax.get_ylim()) )
    min_e  = min( min(ax.get_xlim()) ,  min(ax.get_ylim()) )

    ax.plot([min_e, max_e], [min_e, max_e], color = 'k', linestyle = '--', zorder = 0)

    ax.set_xlabel('DFT Energy (eV/atom)')
    ax.set_ylabel('MLIP Energy (eV/atom)')

    ax.legend(fontsize = 8 , loc = 'upper left',  borderpad = 0.1, handletextpad = 0.3)
    ax.minorticks_on()



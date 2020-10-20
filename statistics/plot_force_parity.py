import numpy as np
from . import struct_colors, dyn_markers, dyn_types, struct_types
from . import read_evaluation_data, get_force_list, collapse_sub_lists


def plot_force_parity(ax, data_sets, 
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
                
                image_pairs = read_evaluation_data(filename = fname, struct_types = [struct_type], dyn_types = [dyn_type])
                if len(image_pairs)>0:
                    cache_forces, data_forces = get_force_list(image_pairs)
                    x, y = collapse_sub_lists(data_forces), collapse_sub_lists(cache_forces)
                    ax.plot(x,y,
                        color = color*lightness,
                        marker = marker,  markersize = 1, linestyle = '',
                        label = label, zorder = zorder)








    ax.set_xlabel('DFT Force (eV/Å)')
    ax.set_ylabel('MLIP Force (eV/Å)')
    
    max_e  = max( max(ax.get_xlim()) ,  max(ax.get_ylim()) )
    min_e  = min( min(ax.get_xlim()) ,  min(ax.get_ylim()) )

    ax.plot([min_e, max_e], [min_e, max_e], color = 'k', linestyle = '--', zorder = 0)
    ax.legend(fontsize = 8 , loc = 'upper left',  borderpad = 0.1)
    ax.minorticks_on()

import numpy as np
from . import struct_colors, dyn_markers, dyn_types, struct_types
from . import read_evaluation_data, get_force_list, compute_force_error_list
from . import compute_rms_force_error_by_atom, compute_rms_force_error_by_image
from . import compute_force_norms_by_image, collapse_sub_lists
from . import compute_force_cosines_by_atom, compute_force_norms_by_atom
from . import nice_bins_percentile
from matplotlib import cm

formula_angle  = r"$\theta_\mathbf{F} = \cos^{-1} \left (   \frac{\mathbf{F}_{MLIP} \cdot \mathbf{F}_{DFT} }{\left | \mathbf{F}_{MLIP}  \right | \left | \mathbf{F}_{DFT}  \right |} \right )$"


def plot_force_angle_heatmaps(axes, 
                data_sets, 
                struct_types = struct_types,
                dyn_types = dyn_types,
                struct_colors = struct_colors,
                dyn_markers = dyn_markers,
                bad_data_traj_list = [],
                cmap = cm.get_cmap('plasma')):
    


    deg = 180/np.pi
    bin_size = 2
    my_bins = np.arange(0,180+bin_size/2, bin_size)

    ylims = (0,180)
    
    from matplotlib.ticker import MultipleLocator
    
    import copy
    cmap_tweaked = copy.copy(cmap)
    from matplotlib.colors import Colormap
    Colormap.set_under(cmap_tweaked, color=(1,1,1,0))
    Colormap.set_over(cmap_tweaked, color=(1,1,1,0))
    
    for di, data_set in enumerate(data_sets):


        fname =  data_set[0]
        data_name = data_set[1]
        #color = np.array(get_color(data_set[2]))
        #lightness = data_set[3]
        #zorder = zorders[di]
        
        image_pairs = read_evaluation_data(filename = fname, struct_types = struct_types, dyn_types = dyn_types)
        cache_forces, data_forces = get_force_list(image_pairs)
        force_cosines_by_atom = compute_force_cosines_by_atom(cache_forces, data_forces)
        data_force_norms_by_atom = compute_force_norms_by_atom(data_forces)

        X = collapse_sub_lists(data_force_norms_by_atom )
        Y = deg*np.arccos( collapse_sub_lists(force_cosines_by_atom))
        xbins = np.linspace(0,X.max(), 100)
        axes[di].hist2d(X, Y, bins = (xbins, my_bins), vmin=1, cmap = cmap_tweaked)
        

        #label = data_name + '\nMean: %.2f°\nRMS: %.2f°'%(mean_force_angles, rms_force_angles)
        
        axes[di].set_title(data_name, fontsize= 8)

        axes[di].set_ylim(ylims )
        #ax2.legend(fontsize = 8, handletextpad = 0.3, borderpad = 0.1)
        axes[di].minorticks_on()
    
        
        axes[di].yaxis.set_major_locator(MultipleLocator(base=30))
        axes[di].yaxis.set_minor_locator(MultipleLocator(base=5))
    
    
        if di == 0:
            axes[di].set_ylabel('Force Angle, '+formula_angle +' (°)')
            axes[di].set_xlabel('DFT Force (eV/Å)')
            #axes[di].set_ylabel('MLIP Force Error (eV/Å)')
    

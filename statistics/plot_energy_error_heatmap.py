import numpy as np
from . import struct_colors, dyn_markers, dyn_types, struct_types
from . import read_evaluation_data
from . import get_energy_lists
#from . import get_force_list, compute_force_error_list
#from . import compute_rms_force_error_by_atom, compute_rms_force_error_by_image
#from . import compute_force_norms_by_image, collapse_sub_lists
#from . import compute_force_cosines_by_atom, compute_force_norms_by_atom
from . import nice_bins_percentile
from matplotlib import cm

formula_angle  = r"$\theta_\mathbf{F} = \cos^{-1} \left (   \frac{\mathbf{F}_{MLIP} \cdot \mathbf{F}_{DFT} }{\left | \mathbf{F}_{MLIP}  \right | \left | \mathbf{F}_{DFT}  \right |} \right )$"


def plot_energy_error_heatmap(ax, 
                data_set, 
                struct_types = struct_types,
                dyn_types = dyn_types,
                bad_data_traj_list = [],
                xbin_size = None,
                ybin_size = None,
                use_meV_y = False,
                use_meV_x = False,
                cmap = cm.get_cmap('plasma')):
    


    scalex = 1
    scaley = 1
    if use_meV_x: scalex = 1000
    if use_meV_y: scaley = 1000
    
    #from matplotlib.ticker import MultipleLocator
    #ax.yaxis.set_major_locator(MultipleLocator(base=30))
    #ax.yaxis.set_minor_locator(MultipleLocator(base=5))
    
    
    import copy
    cmap_tweaked = copy.copy(cmap)
    from matplotlib.colors import Colormap
    Colormap.set_under(cmap_tweaked, color=(1,1,1,0))
    Colormap.set_over(cmap_tweaked, color=(1,1,1,0))
    

    fname =  data_set[0]
    data_name = data_set[1]

    
    image_pairs = read_evaluation_data(filename = fname,
        struct_types = struct_types,
        dyn_types = dyn_types,
        bad_data_traj_list = bad_data_traj_list)
    
    cache_energy, data_energy = get_energy_lists(image_pairs)
    energy_error = cache_energy - data_energy
    error_rmse = np.sqrt(np.mean( energy_error**2)) * scaley
    error_mae  = np.mean(np.absolute(energy_error)) * scaley

    X = data_energy * scalex
    Y = energy_error * scaley
    if xbin_size is None:
        xbins = np.linspace(X.min(),X.max(), 100)
    else:
        xbins = np.arange(X.min(),  X.max(), xbin_size)
    
    if ybin_size is None:
        ybins = np.linspace(Y.min(),Y.max(), 100)
    else:
        ybins = np.arange(Y.min(), Y.max(), ybin_size)


    
    ax.hist2d(X, Y, bins = (xbins, ybins), vmin=1, cmap = cmap_tweaked)
    

    #label = data_name + '\nMean: %.2f°\nRMS: %.2f°'%(mean_force_angles, rms_force_angles)
    
    ax.set_title(data_name , fontsize= 8)

    #scale = 1.05
    #ax.set_ylim(Y.min()*scale, Y.max()*scale )
    #ax.set_xlim(X.min()*scale, X.max()*scale )
    
    #ax2.legend(fontsize = 8, handletextpad = 0.3, borderpad = 0.1)
    ax.minorticks_on()



    return error_rmse, error_mae

    

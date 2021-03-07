import numpy as np
from . import struct_colors, dyn_markers, dyn_types, struct_types
from . import read_evaluation_data, get_force_list, compute_force_error_list
from . import compute_rms_force_error_by_atom, compute_rms_force_error_by_image
from . import compute_force_norms_by_image, collapse_sub_lists
from . import nice_bins_percentile

def plot_force_error_heatmap(ax,
                data_set, 
                struct_types = struct_types,
                dyn_types = dyn_types,
                bad_data_traj_list = [],
                xbin_size = None,
                ybin_size = None,
                use_meV_y = False,
                use_meV_x = False,
                cmap = cm.get_cmap('viridis')):
    
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



    image_pairs = read_evaluation_data(filename = fname, struct_types = struct_types, dyn_types = dyn_types)
    cache_forces, data_forces = get_force_list(image_pairs)
    force_error_list =  compute_force_error_list(cache_forces, data_forces)
    
    rms_force_error_by_atom = compute_rms_force_error_by_atom(cache_forces, data_forces)
    net_rms_force_error_by_atom =  np.sqrt(np.mean( collapse_sub_lists(rms_force_error_by_atom)**2))
    
    rms_force_error_by_image = compute_rms_force_error_by_image(cache_forces, data_forces)/np.sqrt(3)
    net_rms_force_error_by_image = np.sqrt(np.mean( (rms_force_error_by_image)**2))
    image_force_norms = compute_force_norms_by_image(data_forces)/np.sqrt(3)



    X = scalex * image_force_norms  
    Y = scaley * rms_force_error_by_image 
    if xbin_size is None:
        xbins = np.linspace(X.min(),X.max(), 100)
    else:
        xbins = np.arange(X.min(),  X.max(), xbin_size)
    
    if ybin_size is None:
        ybins = np.linspace(Y.min(),Y.max(), 100)
    else:
        ybins = np.arange(Y.min(), Y.max(), ybin_size)


    
    ax.hist2d(X, Y, bins = (xbins, ybins), vmin=1, cmap = cmap_tweaked)
    


    
    ax.set_title(data_name , fontsize= 8)

    ax.minorticks_on()



    return error_rmse, error_mae

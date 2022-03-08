import numpy as np
from . import struct_colors, dyn_markers, dyn_types, struct_types
from . import read_evaluation_data, get_force_list, compute_force_error_list
from . import compute_rms_force_error_by_atom, compute_rms_force_error_by_image
from . import compute_force_norms_by_image, collapse_sub_lists
from . import nice_bins_percentile
from . import compute_force_norms_by_atom
from matplotlib import cm

def plot_force_error_heatmap(ax,
                image_pairs,
                xbin_size = None,
                ybin_size = None,
                use_meV_y = False,
                use_meV_x = False,
                by_atom = False,
                cmap = cm.get_cmap('viridis')):
    
    scalex = 1
    scaley = 1
    if use_meV_x: scalex = 1000
    if use_meV_y: scaley = 1000
    
    
    def pick_bins(abin_size, A):
        if abin_size is None:
            abins = np.linspace(A.min(), A.max(), 100)
        else:
            abins = np.arange(A.min(),  A.max(), abin_size)
        return abins
    #from matplotlib.ticker import MultipleLocator
    #ax.yaxis.set_major_locator(MultipleLocator(base=30))
    #ax.yaxis.set_minor_locator(MultipleLocator(base=5))
    
    
    import copy
    cmap_tweaked = copy.copy(cmap)
    from matplotlib.colors import Colormap
    Colormap.set_under(cmap_tweaked, color=(1,1,1,0))
    Colormap.set_over(cmap_tweaked, color=(1,1,1,0))
    




    #image_pairs = read_evaluation_data(filename = fname, struct_types = struct_types, dyn_types = dyn_types)
    cache_forces, data_forces = get_force_list(image_pairs)
    force_error_list =  compute_force_error_list(cache_forces, data_forces)
    
    if by_atom:
        atom_force_norms = compute_atom_force_norms(data_forces)
        rms_force_error_by_atom = compute_rms_force_error_by_atom(cache_forces, data_forces)
        net_rms_force_error_by_atom =  np.sqrt(np.mean( collapse_sub_lists(rms_force_error_by_atom)**2))
        
        X          = scalex * atom_force_norms
        Y          = scaley * rms_force_error_by_atom
        error_rmse = scaley * net_rms_force_error_by_atom
        
    else:
        image_force_norms = compute_force_norms_by_image(data_forces)/np.sqrt(3)
        rms_force_error_by_image = compute_rms_force_error_by_image(cache_forces, data_forces)/np.sqrt(3)
        net_rms_force_error_by_image = np.sqrt(np.mean( (rms_force_error_by_image)**2))

        X          = scalex * image_force_norms
        Y          = scaley * rms_force_error_by_image
        error_rmse = scaley * net_rms_force_error_by_image
    

    xbins = pick_bins(xbin_size, X)
    ybins = pick_bins(ybin_size, Y)

    ax.hist2d(X, Y, bins = (xbins, ybins), vmin=1, cmap = cmap_tweaked)

    #ax.set_title(data_name , fontsize= 8)

    ax.minorticks_on()

    return error_rmse

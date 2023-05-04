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

#formula_angle  = r"$\theta_\mathbf{F} = \cos^{-1} \left (   \frac{\mathbf{F}_{MLIP} \cdot \mathbf{F}_{DFT} }{\left | \mathbf{F}_{MLIP}  \right | \left | \mathbf{F}_{DFT}  \right |} \right )$"


def plot_energy_error_heatmap(ax, 
                image_pairs,
                xbin_size = None,
                ybin_size = None,
                use_meV_y = False,
                use_meV_x = False,
                cmap = cm.get_cmap('plasma'),
                reference_energies = None):
    '''Plots energy error vs energy heatmap'''


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
    
    if reference_energies is None:
        ref_energies = np.zeros(len(image_pairs))
    else:
        ref_e_per_atom = [ 
            reference_energies[i]/len(image_pairs[i][0]) 
            for i in range(len(image_pairs)) ]
            
    
    
    import copy
    cmap_tweaked = copy.copy(cmap)
    from matplotlib.colors import Colormap
    Colormap.set_under(cmap_tweaked, color=(1,1,1,0))
    Colormap.set_over(cmap_tweaked, color=(1,1,1,0))
    

    
    cache_energy, data_energy = get_energy_lists(image_pairs)
    energy_error = cache_energy - data_energy
    error_rmse = np.sqrt(np.mean( energy_error**2)) * scaley
    error_mae  = np.mean(np.absolute(energy_error)) * scaley

    X = (data_energy - ref_e_per_atom) * scalex
    Y = energy_error * scaley

    
    xbins = pick_bins(xbin_size, X)
    ybins = pick_bins(ybin_size, Y)


    
    ax.hist2d(X, Y, bins = (xbins, ybins), vmin=1, cmap = cmap_tweaked)
    
    #ax.set_title(data_name , fontsize= 8)
    ax.minorticks_on()



    return error_rmse, error_mae

    

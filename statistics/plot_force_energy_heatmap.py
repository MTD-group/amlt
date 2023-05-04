import numpy as np
from . import struct_colors, dyn_markers, dyn_types, struct_types
from . import read_evaluation_data, get_force_list, compute_force_error_list
from . import compute_rms_force_error_by_atom, compute_rms_force_error_by_image
from . import get_energy_lists
from . import compute_force_norms_by_image, collapse_sub_lists
from . import nice_bins_percentile
from . import compute_force_norms_by_atom
from matplotlib import cm

def plot_force_energy_heatmap(ax,
                image_pairs,
                xbin_size = None,
                ybin_size = None,
                use_meV_y = False,
                use_meV_x = False,
                by_atom = False,
                use_logy = True,
                cmap = cm.get_cmap('viridis'),
                reference_energies = None):
    '''Plots forces vs. energy heatmap'''
    
    scalex = 1
    scaley = 1
    if use_meV_x: scalex = 1000
    if use_meV_y: scaley = 1000
    
    
    def pick_bins(abin_size, A, use_log=False):
        if abin_size is None:
            if not use_log:
                abins = np.linspace(A.min(), A.max(), 100)
            else:
                abins = np.logspace(np.log10(A.min()), np.log10(A.max()), 100, base=10)
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
    energies_per_atom = data_energy - ref_e_per_atom

    cache_forces, data_forces = get_force_list(image_pairs)
    
    if by_atom:
        energies_per_atom_expanded = []
        for i in range(len(image_pairs)):
            energies_per_atom_expanded+= len(image_pairs[i][0])*[energies_per_atom[i]]
            
        atom_force_norms = collapse_sub_lists( compute_force_norms_by_atom(data_forces))
        X          = scalex * energies_per_atom_expanded
        Y          = scaley * atom_force_norm
    else:
        image_force_norms = compute_force_norms_by_image(data_forces)/np.sqrt(3)
        X          = scalex * energies_per_atom
        Y          = scaley * image_force_norms

    

    xbins = pick_bins(xbin_size, X)
    ybins = pick_bins(ybin_size, Y, use_log=use_logy)
    print(ybins)
    ax.hist2d(X, Y, bins = (xbins, ybins), vmin=1, cmap = cmap_tweaked)

    if use_logy:
        ax.set_yscale('log')


    ax.minorticks_on()

    #return error_rmse

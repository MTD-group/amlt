import numpy as np
from . import struct_colors, dyn_markers, dyn_types, struct_types
from . import read_evaluation_data, get_force_list, compute_force_error_list
from . import compute_rms_force_error_by_atom, compute_rms_force_error_by_image
from . import compute_force_norms_by_image, collapse_sub_lists
from . import nice_bins_percentile
from . import compute_force_norms_by_atom
from . import compute_force_cosines_by_atom 

from matplotlib import cm

def plot_force_angle_polar_heatmap(ax,
                data_set, 
                struct_types = struct_types,
                dyn_types = dyn_types,
                bad_data_traj_list = [],
                theta_zero_location = "S",
                xbin_size = None,
                ybin_size = None,
                use_meV_y = False,
                #use_meV_x = False,
                by_atom = True,
                cmap = cm.get_cmap('plasma')):
    
    scalex = 1
    scaley = 1
    #if use_meV_x: scalex = 1000
    if use_meV_y: scaley = 1000
    
    deg = 180/np.pi
    
    def pick_bins(abin_size, A):
        if abin_size is None:
            abins = np.linspace(np.min(A), np.max(A), 100)
        else:
            abins = np.arange(np.min(A), np.max(A), abin_size)
        return abins
        
    if xbin_size is None:
            xbin_size = 2
    xbins = np.arange(0,  180, xbin_size)
    
    from matplotlib.ticker import MultipleLocator
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
    
    if by_atom:
        force_norms_by_atom = compute_force_norms_by_atom(data_forces)
        rms_force_error_by_atom = compute_rms_force_error_by_atom(cache_forces, data_forces)
        #net_rms_force_error_by_atom =  np.sqrt(np.mean( collapse_sub_lists(rms_force_error_by_atom)**2))
        force_cosines_by_atom = compute_force_cosines_by_atom(cache_forces, data_forces)
        
        
        force_angles_by_atom = deg*np.arccos(collapse_sub_lists(force_cosines_by_atom))
        mean_force_angles = np.mean(force_angles_by_atom)
        rms_force_angles = np.sqrt(np.mean( (force_angles_by_atom)**2))
        
        X          = scalex * force_angles_by_atom
        Y          = scaley * collapse_sub_lists(rms_force_error_by_atom)
        error_rmse = scaley * rms_force_angles
        error_mae  = scaley * mean_force_angles
        
    else:
        image_force_norms = compute_force_norms_by_image(data_forces)/np.sqrt(3)
        rms_force_error_by_image = compute_rms_force_error_by_image(cache_forces, data_forces)/np.sqrt(3)
        net_rms_force_error_by_image = np.sqrt(np.mean( (rms_force_error_by_image)**2))

        X          = scalex * image_force_norms
        Y          = scaley * rms_force_error_by_image
        error_rmse = scaley * net_rms_force_error_by_image
    

    #xbins = pick_bins(xbin_size, X)
    ybins = pick_bins(ybin_size, Y)

    ax.hist2d(X/deg, Y, bins = (xbins/deg, ybins), vmin=1, cmap = cmap_tweaked)

    ax.set_title(data_name , fontsize= 8)


    #ax_force_polar.legend(fontsize = 8)
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    ax.set_theta_zero_location(theta_zero_location)

    ax.xaxis.set_major_locator(MultipleLocator(base=45/deg))
    #ax.xaxis.set_minor_locator(MultipleLocator(base=15/deg))
    
    dtheta = 15
    theta_grid = np.arange(0,180+dtheta/2, dtheta )
    for theta in theta_grid:
        ax.axvline(theta/deg, color = 'grey', lw = 0.8, zorder = -1)
    #ax.set_thetagrids()

    #ax.minorticks_on()

    return error_rmse, error_mae

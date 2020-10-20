import numpy as np
from . import struct_colors, dyn_markers, dyn_types, struct_types
from . import read_evaluation_data, get_force_list, compute_force_error_list
from . import compute_rms_force_error_by_atom, compute_rms_force_error_by_image
from . import compute_force_norms_by_image, collapse_sub_lists
from . import compute_force_cosines_by_atom, compute_force_norms_by_atom
from . import nice_bins_percentile

formula_angle  = r"$\theta_\mathbf{F} = \cos^{-1} \left (   \frac{\mathbf{F}_{MLIP} \cdot \mathbf{F}_{DFT} }{\left | \mathbf{F}_{MLIP}  \right | \left | \mathbf{F}_{DFT}  \right |} \right )$"

def plot_force_angles(ax1, ax2, 
                data_sets, 
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


    deg = 180/np.pi
    bin_size = 2
    my_bins = np.arange(0,180+bin_size/2, bin_size)


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
                    cache_forces, data_forces = get_force_list(image_pairs)
                    force_cosines_by_atom = compute_force_cosines_by_atom(cache_forces, data_forces)
                    data_force_norms_by_atom = compute_force_norms_by_atom(data_forces)
                    x, y = collapse_sub_lists(data_force_norms_by_atom ), deg*np.arccos( collapse_sub_lists(force_cosines_by_atom)), 
                    ax1.plot(x,y,
                        color = color*lightness,
                        marker = marker,  markersize = 3, linestyle = '',
                        label = label, zorder = zorder)



    #axes2f[0].axhline(0, color = 'grey', linestyle = '--')
    ax1.set_xlabel('DFT Force by Atom (eV/Ang)')
    ax1.set_ylabel('Force Angle, '+formula_angle +' (°)')
    ax1.legend(fontsize = 8,  borderpad = 0.1)
    
    #ylims = auto_limits(test_error)
    #ylims = (-0.3, 0.7)
    #axes[0].set_ylim( ylims)

    ax1.minorticks_on()
    
    #ymax = axes2f[0].get_ylim()[1]
    ax1.set_ylim(0,180)
    ylims = ax1.get_ylim()
    from matplotlib.ticker import MultipleLocator
    ax1.yaxis.set_major_locator(MultipleLocator(base=30))
    ax1.yaxis.set_minor_locator(MultipleLocator(base=5))
    
    xmax = ax1.get_xlim()[1]
    ax1.set_xlim(0,xmax)
    
    
    
    ###### histogram part
    
    for di, data_set in enumerate(data_sets):


        fname =  data_set[0]
        data_name = data_set[1]
        color = np.array(get_color(data_set[2]))
        lightness = data_set[3]
        zorder = zorders[di]
        
        image_pairs = read_evaluation_data(filename = fname, struct_types = struct_types, dyn_types = dyn_types)
        cache_forces, data_forces = get_force_list(image_pairs)
        force_cosines_by_atom = compute_force_cosines_by_atom(cache_forces, data_forces)

        force_angles = deg * np.arccos(collapse_sub_lists(force_cosines_by_atom))
        mean_force_angles = np.mean(force_angles) 
        rms_force_angles  = np.sqrt(np.mean( (force_angles)**2))
        
        #my_bins =  nice_bins_percentile(rms_force_error_by_image, bins_in_window = 30, window_percent = 90)
        label = data_name + '\nMean: %.2f°\nRMS: %.2f°'%(mean_force_angles, rms_force_angles)
        
        ax2.hist(force_angles, my_bins, density = False,
            orientation="horizontal", color = color*lightness, label = label, zorder = zorder )

    ax2.set_ylim(ylims )
    ax2.legend(fontsize = 8, handletextpad = 0.3, borderpad = 0.1)
    ax2.minorticks_on()
    ax2.set_xlabel('# of Samples')
    #pyplot.setp(ax2.get_yticklabels(), visible=False)
    #ax2.yaxis.set_ticklabels([])

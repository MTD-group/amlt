from .analysis import cancel_net_force, read_evaluation_data, get_energy_lists
from .analysis import compute_force_norms_by_atom, get_force_list
from .analysis import compute_force_error_list
from .analysis import get_composition_and_atom_count_arrays
from .analysis import compute_force_norms_by_image, compute_relative_force_error_list
from .analysis import compute_rms_relative_force_error_by_atom
from .analysis import compute_rms_relative_force_error, compute_force_cosines_by_atom
from .analysis import compute_rms_force_error_by_atom
from .analysis import compute_rms_force_error_by_image, collapse_sub_lists
from .analysis import compute_force_cosines_by_atom, compute_force_cosines_by_image



import numpy as np


struct_colors = {
    'random': 'orange',
    'polymorphD3':'teal',
    'known':'purple'}

dyn_markers = {
    'md'   : 'o', # circle
    'relax': 'v', # downward triangle
    'sp'   : '*', # star
    'dimer': '_', # horizontal line
    'ce'   : 'p', # pentagon
    } 
    
    
struct_types = [ 'random', 'polymorphD3', 'known'  ]
dyn_types    = [ 'md','relax','sp','ce','dimer']



def nice_bins_percentile(data, bins_in_window = 30, window_percent = 90):

    window = np.percentile(data, [50 - window_percent/2.0, 50 + window_percent/2.0])
    bin_step = (window.max()-window.min())/bins_in_window
    #rms = sqrt(mean(data**2))
    #bin_step = bin_step_scale*rms
    bin_min = bin_step * np.floor(data.min()/bin_step )
    bin_max = data.max()+bin_step

    bins = np.arange(bin_min,bin_max , bin_step  )
    if len(bins)< 10:
        bins = np.linspace(data.min(), data.max(), 10)
    #print(bins[-1]/bin_step, data.max()/bin_step )
    return bins
    
    
    

from .plot_energy_parity import plot_energy_parity
from .plot_force_parity import plot_force_parity
from .plot_energy_error import plot_energy_error
from .plot_force_error import plot_force_error
from .plot_force_angles import plot_force_angles
from .plot_force_angle_heatmaps import plot_force_angle_heatmaps
from .plot_energy_error_heatmap import plot_energy_error_heatmap

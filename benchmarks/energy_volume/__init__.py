import numpy as np
import os




def read_energy_volume_traj(file_name):
        '''Returns [V/atom,E/atom] sorted by volume/atom from traj file'''
        from ase import io
        traj = io.Trajectory(os.path.abspath(file_name),'r')
        #data = [ (im.get_volume()/len(im), im.get_potential_energy(force_consistent = True)/len(im)) for im in traj]
        data = [ (im.get_volume()/len(im), im.get_potential_energy()/len(im)) for im in traj]
        traj.close()
        data = np.asarray(data).T
        smap = np.argsort(data[0])

        return np.array([data[0][smap],data[1][smap]] )




bulk_modulus_conversion_factor =  (1.6021766208e-19)/(1e-30 * 1e9) # eV to GPa

def format_Murnaghan(name, poptn):
    E0 = poptn[0]
    V0 = poptn[1] # should be equilibrium
    K0 = poptn[2]
    K0_prime = poptn[3]
    nfit = poptn[4]
    entries = (name.ljust(15), 'Equilibrium Volume (Ang^3/atom) %.3f'%V0, 
        'Equilibrium bulk modulus (GPa) %.3f'%(K0*bulk_modulus_conversion_factor),
        'Min Energy %.3f'%E0 , 'Points fit %i'%nfit)
    fit_string = '%s %s %s %s %s'%entries
    print(fit_string)
    #print (name.ljust(15), 'Equilibrium Volume (Ang^3/atom) %.3f'%V0, 'Equilibrium bulk modulus (GPa) %.3f'%(K0*bulk_modulus_conversion_factor), 'Min Energy %.3f'%E0 , 'Points fit %i'%nfit) 

    table_data = [V0, K0*bulk_modulus_conversion_factor, E0*1000 ] # meV
    return table_data, fit_string


def Murnaghan_E_of_V(V, E0, V0, K0, K0_prime):
    E = E0 + K0*V0*( 1.0/(K0_prime*(K0_prime-1.0))*(V/V0)**(1.0 - K0_prime) + (1.0/K0_prime) * (V/V0) - 1.0/(K0_prime-1.0) )
    return E

def fit_Murnaghan(data, K0_guess = 200.0, K0_prime_guess = 3.0,
    trim_high_energy = True, energy_fit_window = 0.2 , V_window = None,
    verbose = False, name = None):
    '''data is [V,E], K0_guess is in GPa,  V_window is the range to look for energy minimum in, and will restrict fitting to this range if high energy is trimmed '''
    V, E = data[0], data[1]
    from scipy.optimize import curve_fit

    if V_window == None:
        V_window = ( V.min(), V.max())
    
    #print(V_window)

    def find_min_E_index(E):
        index_min = 0
        
        for index in range(len(E)):
            if V_window[0] <= V[index] and V[index] <= V_window[1]:
                index_min = index # find soemthing in the correct range

        for index in range(len(E)):
            if E[index] < E[index_min]:
                if V_window[0] <= V[index] and V[index] <= V_window[1]:
                    index_min = index
        return index_min
    

    min_index = find_min_E_index(E)
        
    E_data_min = E[min_index]
    V_data_min = V[min_index]

    
    p0 = [E_data_min, V_data_min, K0_guess/bulk_modulus_conversion_factor, K0_prime_guess ] 
    #print('p0 guess', p0)
    
    if trim_high_energy:
        V_data = []
        E_data = []
        for i in range(len(E)):
            if E[i] < (E_data_min + energy_fit_window):
                if V_window[0] <= V[i] and V[i] <= V_window[1]:
                    E_data.append(E[i])
                    V_data.append(V[i])
    else:
        V_data = V
        E_data = E

    #print('data for fit',V_data, E_data)

    if False:
        popt, pcov = curve_fit(f = Murnaghan_E_of_V, xdata = V_data, ydata = E_data, p0 = p0)
    else:
        try:
            popt, pcov = curve_fit(f = Murnaghan_E_of_V, xdata = V_data, ydata = E_data, p0 = p0)
        except:
            popt = [0,0,0,0]

    poptn = list(popt)+[len(V_data)]
    
    if verbose:
        if name is None:
            mat_name = ''
        else:
            mat_name = name
        format_Murnaghan(name, poptn)
    return poptn



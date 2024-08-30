# KG 7/21/2021

import numpy as np
from code3_PBC_correct import PBC_correct

# A few unit conversions
AMU = 1.6605402e-27  # [kg]
Angstrom = 1.0e-10   # [m]
THz = 1.0e12 # [Hz]
kB = 1.38064852e-23 # J/K

def modal_bottleneck_area(x, L, coeff, w, v, Natoms, mass, NA, NB, T):
    Nmodes = 3 * Natoms
    # initiate a number of variables
    dx_cartesian = np.zeros((Natoms,3))
    bn_area_modal = np.zeros(Nmodes)
    dist_direct = np.zeros([3,3])  # the first 3 corresponds to the number of atoms
    ##### these are for calculating the triangular area
    Natoms_trainagle_calc = 3
    d = np.zeros(Natoms_trainagle_calc)
    #####
    for n in np.arange(0, Nmodes, 1): # looping over all the modes
        w_rad = 2*np.pi*THz*w[n]
        Q = 1./(w_rad) * np.sqrt(kB*T) # formula from Martin Dove lattice dynamics book & https://doi.org/10.1063/1.5081722
        for i in np.arange(0,Natoms,1): # looping over all the atoms
            e_vector = v[3*i+0:3*i+3, n]
            dx_cartesian[i] = 1./np.sqrt(mass[i]*AMU) * Q * e_vector / Angstrom # Angstrom # formula from https://doi.org/10.1063/1.5081722
        # Convert dx from Cartesian to direct
        dx_direct = dx_cartesian @ np.linalg.inv(L)
        x_plus_dx_direct = x + dx_direct
        #### for the displaced structure, calculate the traingular area formed between the the edges formed between 1 nn B-site and 2 nn A-sites
        # first, calculate all distances between the  B-site and 2 A-sites and apply PBC on them
        dist_direct[0] = x_plus_dx_direct[NA[0]] - x_plus_dx_direct[NA[1]]
        dist_direct[1] = x_plus_dx_direct[NA[0]] - x_plus_dx_direct[NB[0]]
        dist_direct[2] = x_plus_dx_direct[NA[1]] - x_plus_dx_direct[NB[0]]
        dist_correcting = PBC_correct(dist_direct, L, coeff, Natoms_trainagle_calc) # apply PBC and convert to cartesian
        dist_correcting.perform_PBC_correct()
        dist_cartesian = dist_correcting.dx_cartesian
        # calclate the triangular area
        count = 0
        for i in np.arange(0,Natoms_trainagle_calc,1):
            d[count] = np.sqrt(np.sum(np.square(dist_cartesian[i])))
            count += 1
        # second, calculate the area using the formula based on having 3 sides 
        # https://www.mathopenref.com/heronsformula.html
        p = np.sum(d)/2
        bn_area_modal[n] = np.sqrt(p*(p-d[0])*(p-d[1])*(p-d[2])) # A^2

    return bn_area_modal

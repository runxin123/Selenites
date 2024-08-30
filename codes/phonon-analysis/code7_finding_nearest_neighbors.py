# KG 7/18/2021

import numpy as np
import code7_finding_nearest_neighbors as nn_find
from code3_PBC_correct import PBC_correct

def nearest_neighbor_find(idx_hopping_ion, atomNs, atomsp, x, L, coeff, Natoms):
    dx_direct = np.empty([Natoms, 3])
    for i in np.arange(0,Natoms,1):
        dx_direct[i] = x[i] - x[idx_hopping_ion]
    dx_correcting = PBC_correct(dx_direct, L, coeff, Natoms) # correcting dx for periodic boundary conditions
    dx_correcting.perform_PBC_correct()
    dist_array = np.empty([Natoms])
    count = 0
    for i in np.arange(0,Natoms,1):
        dist_array[count] = np.sqrt(np.sum(np.square(dx_correcting.dx_cartesian[i])))
        count += 1
    args = np.argsort(dist_array)
    # finding 2 nn B-sites: NB1 & NB2
    count = 0
    NB = [0] * 2
    for i in np.arange(0,Natoms,1):
        if atomNs[0] == 8:           # This is for when we have no substituted atoms # Not good --- hard coded!!
            if args[i]>=atomNs[0] and args[i]<atomNs[0]+atomNs[1]:
                NB[count] = args[i]
                count += 1
                if count>=2:
                    break
        else:
            if args[i]>=atomNs[0]+atomNs[1] and args[i]<atomNs[0]+atomNs[1]+atomNs[2]:
                NB[count] = args[i]
                count += 1
                if count>=2:
                    break

    # finding 4 nn A-sites: NA1, NA2, ...
    NA = [0] * 4
    count = 0
    for i in np.arange(0,Natoms,1):
        if atomNs[0] == 8:           # This is for when we have no substituted atoms
            if args[i]<atomNs[0]:
                NA[count] = args[i]
                count += 1
                if count>=4:
                    break
        else:                        # This is for when we have substituted atoms
            if args[i]<atomNs[0]+atomNs[1]:
                NA[count] = args[i]
                count += 1
                if count>=4:
                    break
            
    return NA, NB

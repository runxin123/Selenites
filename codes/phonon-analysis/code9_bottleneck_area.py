# KG 7/20/2021

import numpy as np
from code3_PBC_correct import PBC_correct
import pymatgen.core as pmg

def bottleneck_area_calc(dx, x, L, coeff, Natoms, atomNs, atomsp):
    # find the index of the hopping ion
    maxdx_idx = 0
    maxdx = 0
    for i in np.arange(0,Natoms,1):
        temp = np.sqrt(np.sum(np.square(dx[i])))
        if temp > maxdx:
            maxdx_idx = i
            maxdx = temp
    idx_hopping_ion = maxdx_idx

    # calculate all distances of all atoms from idx_hopping_ion and apply PBC on them
    dist_direct = x - x[idx_hopping_ion]
    dist_correcting = PBC_correct(dist_direct, L, coeff, Natoms)
    dist_correcting.perform_PBC_correct()
    dist_cartesian = dist_correcting.dx_cartesian

    # Get the 1 nn B-site and 2 nn A-sites
    howmany_NBs = 1
    howmany_NAs = 2
    ####---####
    dist_array = np.empty([Natoms])
    count = 0
    for i in np.arange(0,Natoms,1):
        dist_array[count] = np.sqrt(np.sum(np.square(dist_cartesian[i])))
        count += 1
    args = np.argsort(dist_array)

    if len(atomsp) == 3:
        # finding nn B-site
        count = 0
        NB = [0] * howmany_NBs
        for i in np.arange(0,Natoms,1):
            if args[i]>=atomNs[0] and args[i]<atomNs[0]+atomNs[1]:
                NB[count] = args[i]
                count += 1
                if count >= howmany_NBs:
                    break
        # finding nn A-sites
        count = 0
        NA = [0] * howmany_NAs
        for i in np.arange(0,Natoms,1):
            if args[i]<atomNs[0]:
                NA[count] = args[i]
                count += 1
                if count >= howmany_NAs:
                    break
    elif len(atomsp) == 4:
        # finding nn B-site
        count = 0
        NB = [0] * howmany_NBs
        for i in np.arange(0,Natoms,1):
            if args[i]>=atomNs[0]+atomNs[1] and args[i]<atomNs[0]+atomNs[1]+atomNs[2]:
                NB[count] = args[i]
                count += 1
                if count >= howmany_NBs:
                    break
        # finding nn A-sites
        count = 0
        NA = [0] * howmany_NAs
        for i in np.arange(0,Natoms,1):
            if args[i]<atomNs[0]+atomNs[1]:
                NA[count] = args[i]
                count += 1
                if count >= howmany_NAs:
                    break

    #### Calculate the traingular area formed between the the edges formed between 1 nn B-site and 2 nn A-sites
    # first, calculate all distances between the  B-site and 2 A-sites and apply PBC on them
    dist_direct = np.zeros([3,3])  # the first 3 corresponds to the number of atoms
    dist_direct[0] = x[NA[0]] - x[NA[1]]
    dist_direct[1] = x[NA[0]] - x[NB[0]]
    dist_direct[2] = x[NA[1]] - x[NB[0]]
    Natoms_trainagle_calc =3
    dist_correcting = PBC_correct(dist_direct, L, coeff, Natoms_trainagle_calc)
    dist_correcting.perform_PBC_correct()
    dist_cartesian = dist_correcting.dx_cartesian
    d = np.zeros(Natoms_trainagle_calc)
    count = 0
    for i in np.arange(0,Natoms_trainagle_calc,1):
        d[count] = np.sqrt(np.sum(np.square(dist_cartesian[i])))
        count += 1
    # second, calculate the area using the formula based on having 3 sides 
    # https://www.mathopenref.com/heronsformula.html
    p = np.sum(d)/2
    area = np.sqrt(p*(p-d[0])*(p-d[1])*(p-d[2])) # A^2

    return area, NA, NB

def a_few_bond_calc(dx, x1, x3, x4, L, coeff, Natoms, atomNs, atomsp):
    # find the index of the hopping ion
    maxdx_idx = 0
    maxdx = 0
    for i in np.arange(0,Natoms,1):
        temp = np.sqrt(np.sum(np.square(dx[i])))
        if temp > maxdx:
            maxdx_idx = i
            maxdx = temp
    idx_hopping_ion = maxdx_idx

    ###############
    ############### A-O & B-O initial
    ###############
    # calculate all distances of all atoms from idx_hopping_ion and apply PBC on them
    dist_direct = x1 - x1[idx_hopping_ion]
    dist_correcting = PBC_correct(dist_direct, L, coeff, Natoms)
    dist_correcting.perform_PBC_correct()
    dist_cartesian = dist_correcting.dx_cartesian

    # Get the nn B-site and nn A-sites
    howmany_NBs = 2
    howmany_NAs = 4
    ####---####
    dist_array = np.empty([Natoms])
    count = 0
    for i in np.arange(0,Natoms,1):
        dist_array[count] = np.sqrt(np.sum(np.square(dist_cartesian[i])))
        count += 1
    args = np.argsort(dist_array)

    if len(atomsp) == 3:
        # finding nn B-site
        count = 0
        NB = [0] * howmany_NBs
        for i in np.arange(0,Natoms,1):
            if args[i]>=atomNs[0] and args[i]<atomNs[0]+atomNs[1]:
                NB[count] = args[i]
                count += 1
                if count >= howmany_NBs:
                    break
        # finding nn A-sites
        count = 0
        NA = [0] * howmany_NAs
        for i in np.arange(0,Natoms,1):
            if args[i]<atomNs[0]:
                NA[count] = args[i]
                count += 1
                if count >= howmany_NAs:
                    break
    elif len(atomsp) == 4:
        # finding nn B-site
        count = 0
        NB = [0] * howmany_NBs
        for i in np.arange(0,Natoms,1):
            if args[i]>=atomNs[0]+atomNs[1] and args[i]<atomNs[0]+atomNs[1]+atomNs[2]:
                NB[count] = args[i]
                count += 1
                if count >= howmany_NBs:
                    break
        # finding nn A-sites
        count = 0
        NA = [0] * howmany_NAs
        for i in np.arange(0,Natoms,1):
            if args[i]<atomNs[0]+atomNs[1]:
                NA[count] = args[i]
                count += 1
                if count >= howmany_NAs:
                    break

    dist_direct = np.zeros([howmany_NBs,3])  # the first 3 corresponds to the number of atoms
    for i in np.arange(0, howmany_NBs, 1):
        dist_direct[i] = x1[NB[i]] - x1[idx_hopping_ion]
    dist_correcting = PBC_correct(dist_direct, L, coeff, howmany_NBs)
    dist_correcting.perform_PBC_correct()
    dist_cartesian = dist_correcting.dx_cartesian
    d = np.zeros(howmany_NBs)
    for i in np.arange(0, howmany_NBs, 1):
        d[i] = np.sqrt(np.sum(np.square(dist_cartesian[i])))

    B_O2i = np.mean(d)

    dist_direct = np.zeros([howmany_NAs,3])  # the first 3 corresponds to the number of atoms
    for i in np.arange(0, howmany_NAs, 1):
        dist_direct[i] = x1[NA[i]] - x1[idx_hopping_ion]
    dist_correcting = PBC_correct(dist_direct, L, coeff, howmany_NAs)
    dist_correcting.perform_PBC_correct()
    dist_cartesian = dist_correcting.dx_cartesian
    d = np.zeros(howmany_NAs)
    for i in np.arange(0, howmany_NAs, 1):
        d[i] = np.sqrt(np.sum(np.square(dist_cartesian[i])))
    A_O4i = np.mean(d)

    ###############
    ############### A-O & B-O final
    ###############
    # calculate all distances of all atoms from idx_hopping_ion and apply PBC on them
    dist_direct = x4 - x4[idx_hopping_ion]
    dist_correcting = PBC_correct(dist_direct, L, coeff, Natoms)
    dist_correcting.perform_PBC_correct()
    dist_cartesian = dist_correcting.dx_cartesian

    # Get the nn B-site and nn A-sites
    howmany_NBs = 2
    howmany_NAs = 4
    ####---####
    dist_array = np.empty([Natoms])
    count = 0
    for i in np.arange(0,Natoms,1):
        dist_array[count] = np.sqrt(np.sum(np.square(dist_cartesian[i])))
        count += 1
    args = np.argsort(dist_array)

    if len(atomsp) == 3:
        # finding nn B-site
        count = 0
        NB = [0] * howmany_NBs
        for i in np.arange(0,Natoms,1):
            if args[i]>=atomNs[0] and args[i]<atomNs[0]+atomNs[1]:
                NB[count] = args[i]
                count += 1
                if count >= howmany_NBs:
                    break
        # finding nn A-sites
        count = 0
        NA = [0] * howmany_NAs
        for i in np.arange(0,Natoms,1):
            if args[i]<atomNs[0]:
                NA[count] = args[i]
                count += 1
                if count >= howmany_NAs:
                    break
    elif len(atomsp) == 4:
        # finding nn B-site
        count = 0
        NB = [0] * howmany_NBs
        for i in np.arange(0,Natoms,1):
            if args[i]>=atomNs[0]+atomNs[1] and args[i]<atomNs[0]+atomNs[1]+atomNs[2]:
                NB[count] = args[i]
                count += 1
                if count >= howmany_NBs:
                    break
        # finding nn A-sites
        count = 0
        NA = [0] * howmany_NAs
        for i in np.arange(0,Natoms,1):
            if args[i]<atomNs[0]+atomNs[1]:
                NA[count] = args[i]
                count += 1
                if count >= howmany_NAs:
                    break

    dist_direct = np.zeros([howmany_NBs,3])  # the first 3 corresponds to the number of atoms
    for i in np.arange(0, howmany_NBs, 1):
        dist_direct[i] = x4[NB[i]] - x4[idx_hopping_ion]
    dist_correcting = PBC_correct(dist_direct, L, coeff, howmany_NBs)
    dist_correcting.perform_PBC_correct()
    dist_cartesian = dist_correcting.dx_cartesian
    d = np.zeros(howmany_NBs)
    for i in np.arange(0, howmany_NBs, 1):
        d[i] = np.sqrt(np.sum(np.square(dist_cartesian[i])))
    B_O2f = np.mean(d)

    dist_direct = np.zeros([howmany_NAs,3])  # the first 3 corresponds to the number of atoms
    for i in np.arange(0, howmany_NAs, 1):
        dist_direct[i] = x4[NA[i]] - x4[idx_hopping_ion]
    dist_correcting = PBC_correct(dist_direct, L, coeff, howmany_NAs)
    dist_correcting.perform_PBC_correct()
    dist_cartesian = dist_correcting.dx_cartesian
    d = np.zeros(howmany_NAs)
    for i in np.arange(0, howmany_NAs, 1):
        d[i] = np.sqrt(np.sum(np.square(dist_cartesian[i])))
    A_O4f = np.mean(d)

    ###############
    ############### A-O & B-O TS
    ###############
    # calculate all distances of all atoms from idx_hopping_ion and apply PBC on them
    dist_direct = x3 - x3[idx_hopping_ion]
    dist_correcting = PBC_correct(dist_direct, L, coeff, Natoms)
    dist_correcting.perform_PBC_correct()
    dist_cartesian = dist_correcting.dx_cartesian

    # Get the nn B-site and nn A-sites
    howmany_NBs = 1
    howmany_NAs = 1
    ####---####
    dist_array = np.empty([Natoms])
    count = 0
    for i in np.arange(0,Natoms,1):
        dist_array[count] = np.sqrt(np.sum(np.square(dist_cartesian[i])))
        count += 1
    args = np.argsort(dist_array)

    if len(atomsp) == 3:
        # finding nn B-site
        count = 0
        NB = [0] * howmany_NBs
        for i in np.arange(0,Natoms,1):
            if args[i]>=atomNs[0] and args[i]<atomNs[0]+atomNs[1]:
                NB[count] = args[i]
                count += 1
                if count >= howmany_NBs:
                    break
        # finding nn A-sites
        count = 0
        NA = [0] * howmany_NAs
        for i in np.arange(0,Natoms,1):
            if args[i]<atomNs[0]:
                NA[count] = args[i]
                count += 1
                if count >= howmany_NAs:
                    break
    elif len(atomsp) == 4:
        # finding nn B-site
        count = 0
        NB = [0] * howmany_NBs
        for i in np.arange(0,Natoms,1):
            if args[i]>=atomNs[0]+atomNs[1] and args[i]<atomNs[0]+atomNs[1]+atomNs[2]:
                NB[count] = args[i]
                count += 1
                if count >= howmany_NBs:
                    break
        # finding nn A-sites
        count = 0
        NA = [0] * howmany_NAs
        for i in np.arange(0,Natoms,1):
            if args[i]<atomNs[0]+atomNs[1]:
                NA[count] = args[i]
                count += 1
                if count >= howmany_NAs:
                    break

    dist_direct = np.zeros([howmany_NBs,3])  # the first 3 corresponds to the number of atoms
    for i in np.arange(0, howmany_NBs, 1):
        dist_direct[i] = x3[NB[i]] - x3[idx_hopping_ion]
    dist_correcting = PBC_correct(dist_direct, L, coeff, howmany_NBs)
    dist_correcting.perform_PBC_correct()
    dist_cartesian = dist_correcting.dx_cartesian
    d = np.zeros(howmany_NBs)
    for i in np.arange(0, howmany_NBs, 1):
        d[i] = np.sqrt(np.sum(np.square(dist_cartesian[i])))
    B_O_TS = np.mean(d)

    dist_direct = np.zeros([howmany_NAs,3])  # the first 3 corresponds to the number of atoms
    for i in np.arange(0, howmany_NAs, 1):
        dist_direct[i] = x3[NA[i]] - x3[idx_hopping_ion]
    dist_correcting = PBC_correct(dist_direct, L, coeff, howmany_NAs)
    dist_correcting.perform_PBC_correct()
    dist_cartesian = dist_correcting.dx_cartesian
    d = np.zeros(howmany_NAs)
    for i in np.arange(0, howmany_NAs, 1):
        d[i] = np.sqrt(np.sum(np.square(dist_cartesian[i])))
    A_O_TS = np.mean(d)

    atomR = np.zeros(len(atomsp))
    for i in np.arange(0, len(atomsp), 1):
        atomR[i] = float(pmg.Element(atomsp[i]).atomic_radius)

    if len(atomsp) == 3:
        Rcation_avg_weighted_on_stoichiometry = atomR[0]
        Rdopant = None
        RB = atomR[1]
        RA = atomR[0]
    elif len(atomsp) == 4:
        Rcation_avg_weighted_on_stoichiometry = (atomNs[0]*atomR[0] + atomNs[1]*atomR[1])/(atomNs[0]+atomNs[1])
        Rdopant = atomR[1]
        RB = atomR[2]
        RA = atomR[0]
    else:
        print("Something is wrong!!!!!!!")

    return Rcation_avg_weighted_on_stoichiometry, Rdopant, RB, RA, B_O2i, B_O2f, A_O4i, A_O4f, B_O_TS


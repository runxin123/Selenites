# KG 7/18/2021

# formula for Mean-Square-Displacement-matrix is from:
# https://phonopy.github.io/phonopy/formulation.html#thermal-displacement
# After you find the Mean-Square-Displacement-matrix, go and find the eigenvalues of it.
# Those will be the magnitudes of vibration along the prinicipal axes.
# formula for TE magnitude along a speciic axis is obtained from:
# https://phonopy.github.io/phonopy/formulation.html#thermal-displacement

import numpy as np
import code7_finding_nearest_neighbors as nn_find
from math import exp
from scipy import linalg as LA
import yaml
import os
from shutil import copy2

# A few unit conversions
AMU = 1.6605402e-27  # [kg]
Angstrom = 1.0e-10   # [m]
THz = 1.0e12 # [Hz]
hbar = 6.63e-34/(2*np.pi) # J.s
kB = 1.38064852e-23 # J/K

def perform_TE_phonopy(Phonopy_path, fmin, MP):
    # This portion of the code depends on the fact that force constant file has been brought into the folder by the code2 in a previous step.
    os.chdir('./1_phonopy_calc_room')
    #copy2('./tdisp_0.conf','./tdisp.conf') # We will add the hopping pathway direction (coordinates) to the end of tdisp.conf

    with open("./tdisp_0.conf", "r") as f:   ####
        lines = f.readlines()
    for i in range(0, len(lines)):
        templ = lines[i].split()
        if templ[0].lower() == 'fmin':
            lines[i] = 'FMIN = ' + str(fmin) + '\n'
        if templ[0].lower() == 'mp':
            lines[i] = 'MP = ' + MP + '\n'
    filename = "tdisp.conf"
    f_tdisp = open(filename, "w")
    for i in range(0, len(lines)):
        f_tdisp.write("%s\n" %(lines[i][:-1]))
    f_tdisp.close()
    
    os.system(Phonopy_path + "phonopy --readfc tdisp.conf > out.txt")
    os.chdir('..')

    with open("./1_phonopy_calc_room/thermal_displacements.yaml") as f1:
        TE_data_non_proj = yaml.load(f1, Loader=yaml.FullLoader)

    return TE_data_non_proj

def therm_eps_mag_and_projected_hopping_ion(dx, TE_data_non_proj, mass, Natoms, T, fmin, Phonopy_path): ### temperature should be added here
    # find the index of the hopping ion
    maxdx_idx = 0
    maxdx = 0
    for i in np.arange(0,Natoms,1):
        temp = np.sqrt(np.sum(np.square(dx[i])))
        if temp > maxdx:
            maxdx_idx = i
            maxdx = temp
    idx_hopping_ion = maxdx_idx
    print(idx_hopping_ion)
    nhat = dx[idx_hopping_ion] / np.sqrt(np.sum(np.square(dx)))
    os.chdir('./1_phonopy_calc_room')
    file_object = open('tdisp.conf', 'a')
    str2beadded = 'PROJECTION_DIRECTION = %f %f %f' %(nhat[0], nhat[1],nhat[2])
    file_object.write(str2beadded)
    file_object.close()
    os.system(Phonopy_path + "phonopy --readfc tdisp.conf > out.txt")
    os.chdir('..')

    with open("./1_phonopy_calc_room/thermal_displacements.yaml") as f:
        TE_data_proj = yaml.load(f, Loader=yaml.FullLoader)

    '''Nmodes = 3*Natoms
    N_unitcells = 1
    U_cart = np.zeros((3,3), dtype=complex)
    u_proj_squared = 0 # Gamma-point lattice-dynamics calculation
    atom_mass = mass[idx_hopping_ion]*AMU

    
    N_qpoints = len(data['phonon'])
    N_modes_each_qpoints = len(data['phonon'][0]['band']) # Think of this as Nband = Nmodes
    if (N_modes_each_qpoints != Nmodes):
        print('Waaaaaaarnning!!!!!!!!!!!!')
    for nq in np.arange(0, N_qpoints, 1):
        # Let's prepare the freq vector and ev matrix for each q-point first 
        w = np.zeros([N_modes_each_qpoints]) # Initiating the freq array
        v = np.zeros([N_modes_each_qpoints, N_modes_each_qpoints], dtype=complex) # Initiating the ev array
        for nm in np.arange(0, N_modes_each_qpoints, 1): # This loop is not needed. It will be done inside each q-point calculation.
            w[nm] = data['phonon'][nq]['band'][nm]['frequency']
            for na in np.arange(0, Natoms, 1):
                v.real[na*3+0, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][0][0]
                v.real[na*3+1, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][1][0]
                v.real[na*3+2, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][2][0]
                v.imag[na*3+0, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][0][1]
                v.imag[na*3+1, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][1][1]
                v.imag[na*3+2, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][2][1]
        # Then, let's loop over different modes belonging to each q-point
        # nm and n are essentially counters for the # of bands (modes)
        N_qpoints_count = 0
        for n in np.arange(0,Nmodes,1): # looping over all the modes
          w_rad = 2*np.pi*THz*w[n]
          #if (np.absolute(w[n])>fmin):
          if (w[n]>fmin):   # This produces the same result as phonopy
            N_qpoints_count += 1
            nBE = 1./(exp(hbar*w_rad/(kB*T))-1)
            e_vector = v[3*idx_hopping_ion+0:3*idx_hopping_ion+3, n]
            U_cart += hbar/(2*N_unitcells*atom_mass)*1./(w_rad)*(1+2*nBE)*np.outer(np.transpose(e_vector),np.conjugate(np.transpose(e_vector)))
            nhat_dot_ev = np.dot(nhat,e_vector)
            u_proj_squared += hbar/(2*N_unitcells*atom_mass)*1./(w_rad)*(1+2*nBE)*np.square(np.absolute(nhat_dot_ev))

    U_cart /= N_qpoints

    #evalue, evector = np.linalg.eig(U_cart)
    evalue, evector = LA.eigh(U_cart)
    evalue = np.sqrt(np.sort(evalue)) / Angstrom
    hopping_ion_TE_max_evalue = evalue[2]
    
    eff_disp = np.power(evalue[0]*evalue[1]*evalue[2],1./3)
    u_proj = np.sqrt(u_proj_squared) / Angstrom

    r1 = (evalue[1]-evalue[0])/(evalue[2]-evalue[0])
    r2 = evalue[2]/evalue[0]'''

    evalue = np.zeros([3])
    evalue = TE_data_non_proj['thermal_displacements'][0]['displacements'][idx_hopping_ion]

    evalue = np.sqrt(np.sort(evalue))
    print(evalue)
    hopping_ion_TE_max_evalue = evalue[2] 
    r1 = (evalue[1]-evalue[0])/(evalue[2]-evalue[0])
    r2 = evalue[2]/evalue[0]

    eff_disp = np.power(evalue[0]*evalue[1]*evalue[2],1./3)
    u_proj = np.sqrt(TE_data_proj['thermal_displacements'][0]['displacements'][idx_hopping_ion][0])

    return idx_hopping_ion, eff_disp, u_proj, hopping_ion_TE_max_evalue, r1, r2

def therm_eps_mag_nearest_neighbors_to_hopping_ion(idx_hopping_ion, atomNs, atomsp, x, L, coeff, TE_data_non_proj, mass, Natoms, T, fmin):
    # find the 4 nearest A-sites and 2 nearest B-sites to the hopping O atom
    [NA, NB] = nn_find.nearest_neighbor_find(idx_hopping_ion, atomNs, atomsp, x, L, coeff, Natoms)

    '''Nmodes = 3*Natoms
    N_unitcells = 1

    NA_U_cart = np.zeros((len(NA),3,3), dtype=complex)
    NB_U_cart = np.zeros((len(NB),3,3), dtype=complex)'''
    evalue_A = np.zeros((len(NA),3))
    evalue_B = np.zeros((len(NB),3))
    eff_disp_A = np.zeros(len(NA))
    eff_disp_B = np.zeros(len(NB))

    '''N_qpoints = len(data['phonon'])
    N_modes_each_qpoints = len(data['phonon'][0]['band']) # Think of this as Nband = Nmodes
    if (N_modes_each_qpoints != Nmodes):
        print('Waaaaaaarnning!!!!!!!!!!!!')
    for nq in np.arange(0, N_qpoints, 1):
        # Let's prepare the freq vector and ev matrix for each q-point first 
        w = np.zeros([N_modes_each_qpoints]) # Initiating the freq array
        v = np.zeros([N_modes_each_qpoints, N_modes_each_qpoints], dtype=complex) # Initiating the ev array
        for nm in np.arange(0, N_modes_each_qpoints, 1): # This loop is not needed. It will be done inside each q-point calculation.
            w[nm] = data['phonon'][nq]['band'][nm]['frequency']
            for na in np.arange(0, Natoms, 1):
                v.real[na*3+0, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][0][0]
                v.real[na*3+1, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][1][0]
                v.real[na*3+2, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][2][0]
                v.imag[na*3+0, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][0][1]
                v.imag[na*3+1, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][1][1]
                v.imag[na*3+2, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][2][1]
        # Then, let's loop over different modes belonging to each q-point
        # nm and n are essentially counters for the # of bands (modes)
        for n in np.arange(0,Nmodes,1): # looping over all the modes
          w_rad = 2*np.pi*THz*w[n]
          #if (np.absolute(w[n])>fmin):
          if (w[n]>fmin):   # This produces the same result as phonopy
            nBE = 1./(exp(hbar*w_rad/(kB*T))-1)
            for i in np.arange(0,len(NA),1):
                atom_mass = mass[NA[i]]*AMU
                e_vector = v[3*NA[i]+0:3*NA[i]+3, n]
                #if (nq == 2 and n == 10):
                #    print(e_vector)
                NA_U_cart[i] += hbar/(2*N_unitcells*atom_mass)*1./(w_rad)*(1+2*nBE)*np.outer(np.transpose(e_vector),np.conjugate(np.transpose(e_vector)))
            for i in np.arange(0,len(NB),1):
                atom_mass = mass[NB[i]]*AMU
                e_vector = v[3*NB[i]+0:3*NB[i]+3, n]
                NB_U_cart[i] += hbar/(2*N_unitcells*atom_mass)*1./(w_rad)*(1+2*nBE)*np.outer(np.transpose(e_vector),np.conjugate(np.transpose(e_vector)))

    # Normalize it by the number of your discretization (i.e., q-points)
    NA_U_cart /= N_qpoints
    NB_U_cart /= N_qpoints'''

    evalue = np.zeros([3])
    for i in np.arange(0,len(NA),1):
        #evalue_A[i], evector = np.linalg.eig(NA_U_cart[i])
        #evalue_A[i], evector = LA.eigh(NA_U_cart[i])
        evalue_A[i] = TE_data_non_proj['thermal_displacements'][0]['displacements'][NA[i]]
        evalue_A[i] = np.sqrt(np.sort(evalue_A[i]))
        eff_disp_A[i] = np.power(evalue_A[i][0]*evalue_A[i][1]*evalue_A[i][2],1./3)
    for i in np.arange(0,len(NB),1):
        #evalue_B[i], evector = np.linalg.eig(NB_U_cart[i])
        #evalue_B[i], evector = LA.eigh(NB_U_cart[i])
        evalue_B[i] = TE_data_non_proj['thermal_displacements'][0]['displacements'][NB[i]]
        evalue_B[i] = np.sqrt(np.sort(evalue_B[i]))
        eff_disp_B[i] = np.power(evalue_B[i][0]*evalue_B[i][1]*evalue_B[i][2],1./3)

    A_site_TE_mag = np.mean(eff_disp_A)
    B_site_TE_mag = np.mean(eff_disp_B)

    return A_site_TE_mag, B_site_TE_mag

############################### NEW
def therm_eps_mag_average_all_species(atomNs, atomsp, x, L, coeff, TE_data_non_proj, mass, Natoms, T, fmin, AB_first, compound_name, persite_amp_True):
    #define the A, B & O indices -- there is no need to define this.
    #if atomNs[0] == 8: # This is for when we have no substituted atoms
    #    idx_A1 = 0
    #    idx_B1 = 

    if AB_first == 'A' and len(atomNs) == 3:
        endofA = atomNs[0] # A-atoms would only include the main A-site
        endofB = atomNs[0]+atomNs[1]
        endofO = atomNs[0]+atomNs[1]+atomNs[2]
        NA = np.arange(0, endofA, 1)
        NB = np.arange(endofA, endofB, 1)
        NO = np.arange(endofB, endofO, 1)
    elif AB_first == 'A' and len(atomNs) == 4:
        endofA = atomNs[0] +atomNs[1] # A-atoms would include the main A-site and its dopant
        endofB = atomNs[0]+atomNs[1]+atomNs[2]
        endofO = atomNs[0]+atomNs[1]+atomNs[2]+atomNs[3]
        NA = np.arange(0, endofA, 1)
        NB = np.arange(endofA, endofB, 1)
        NO = np.arange(endofB, endofO, 1)
    elif AB_first == 'A' and len(atomNs) == 5:
        endofA = atomNs[0] +atomNs[1] # A-atoms would include the main A-site and its dopant
        endofB = atomNs[0]+atomNs[1]+atomNs[2]+atomNs[3]
        endofO = atomNs[0]+atomNs[1]+atomNs[2]+atomNs[3]+atomNs[4]
        NA = np.arange(0, endofA, 1)
        NB = np.arange(endofA, endofB, 1)
        NO = np.arange(endofB, endofO, 1)
    elif AB_first == 'B' and len(atomNs) == 3:
        endofB = atomNs[0]
        endofA = atomNs[0]+atomNs[1] # A-atoms would only include the main A-site
        endofO = atomNs[0]+atomNs[1]+atomNs[2]
        NB = np.arange(0, endofB, 1)
        NA = np.arange(endofB, endofA, 1)
        NO = np.arange(endofA, endofO, 1)
    elif AB_first == 'B' and len(atomNs) == 4:
        endofB = atomNs[0]+atomNs[1]
        endofA = atomNs[0]+atomNs[1]+atomNs[2] # A-atoms would include the main A-site and its dopant
        endofO = atomNs[0]+atomNs[1]+atomNs[2]+atomNs[3]
        NB = np.arange(0, endofB, 1)
        NA = np.arange(endofB, endofA, 1)
        NO = np.arange(endofA, endofO, 1)

    '''Nmodes = 3*Natoms
    N_unitcells = 1

    # The effect of each phonon will be added to the following matrices
    # Phonons will loop over different q-points and bands (modes)
    NA_U_cart = np.zeros((len(NA),3,3), dtype=complex)
    NB_U_cart = np.zeros((len(NB),3,3), dtype=complex)
    NO_U_cart = np.zeros((len(NO),3,3), dtype=complex)'''

    evalue_A = np.zeros((len(NA),3))
    evalue_B = np.zeros((len(NB),3))
    evalue_O = np.zeros((len(NO),3))
    evalue_all_atoms = np.zeros(3)

    eff_disp_A = np.zeros(len(NA))
    eff_disp_B = np.zeros(len(NB))
    eff_disp_O = np.zeros(len(NO))
    #eff_disp_all_atoms = np.zeros(Natoms)
    #print(Natoms)

    '''N_qpoints = len(data['phonon'])
    N_modes_each_qpoints = len(data['phonon'][0]['band']) # Think of this as Nband = Nmodes
    if (N_modes_each_qpoints != Nmodes):
        print('Waaaaaaarnning!!!!!!!!!!!!')
    for nq in np.arange(0, N_qpoints, 1):
        # Let's prepare the freq vector and ev matrix for each q-point first 
        w = np.zeros([N_modes_each_qpoints]) # Initiating the freq array
        v = np.zeros([N_modes_each_qpoints, N_modes_each_qpoints], dtype=complex) # Initiating the ev array
        for nm in np.arange(0, N_modes_each_qpoints, 1): # This loop is not needed. It will be done inside each q-point calculation.
            w[nm] = data['phonon'][nq]['band'][nm]['frequency']
            for na in np.arange(0, Natoms, 1):
                v.real[na*3+0, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][0][0]
                v.real[na*3+1, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][1][0]
                v.real[na*3+2, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][2][0]
                v.imag[na*3+0, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][0][1]
                v.imag[na*3+1, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][1][1]
                v.imag[na*3+2, nm] = data['phonon'][nq]['band'][nm]['eigenvector'][na][2][1]
        # Then, let's loop over different modes belonging to each q-point
        # nm and n are essentially counters for the # of bands (modes)
        for n in np.arange(0,Nmodes,1): # looping over all the modes
          w_rad = 2*np.pi*THz*w[n]
          #if (np.absolute(w[n])>fmin):
          if (w[n]>fmin):   # This produces the same result as phonopy
            nBE = 1./(exp(hbar*w_rad/(kB*T))-1)
            for i in np.arange(0,len(NA),1):
                atom_mass = mass[NA[i]]*AMU
                e_vector = v[3*NA[i]+0:3*NA[i]+3, n]
                #if (nq == 2 and n == 10):
                #    print(e_vector)
                NA_U_cart[i] += hbar/(2*N_unitcells*atom_mass)*1./(w_rad)*(1+2*nBE)*np.outer(np.transpose(e_vector),np.conjugate(np.transpose(e_vector)))
            for i in np.arange(0,len(NB),1):
                atom_mass = mass[NB[i]]*AMU
                e_vector = v[3*NB[i]+0:3*NB[i]+3, n]
                NB_U_cart[i] += hbar/(2*N_unitcells*atom_mass)*1./(w_rad)*(1+2*nBE)*np.outer(np.transpose(e_vector),np.conjugate(np.transpose(e_vector)))
            for i in np.arange(0,len(NO),1):
                atom_mass = mass[NO[i]]*AMU
                e_vector = v[3*NO[i]+0:3*NO[i]+3, n]
                NO_U_cart[i] += hbar/(2*N_unitcells*atom_mass)*1./(w_rad)*(1+2*nBE)*np.outer(np.transpose(e_vector),np.conjugate(np.transpose(e_vector)))

    # Normalize it by the number of your discretization (i.e., q-points)
    NA_U_cart /= N_qpoints
    NB_U_cart /= N_qpoints
    NO_U_cart /= N_qpoints'''
    np.warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)
    for i in np.arange(0,len(NA),1):
        #evalue_A[i], evector = np.linalg.eig(NA_U_cart[i])
        #evalue_A[i], evector = LA.eigh(NA_U_cart[i])
        #print(i+1, 'A', evalue_A[i][0]/Angstrom**2, evalue_A[i][1]/Angstrom**2, evalue_A[i][2]/Angstrom**2)
        evalue_A[i] = TE_data_non_proj['thermal_displacements'][0]['displacements'][NA[i]]
        evalue_A[i] = np.sqrt(np.sort(evalue_A[i]))
        eff_disp_A[i] = np.power(evalue_A[i][0]*evalue_A[i][1]*evalue_A[i][2],1./3)
    for i in np.arange(0,len(NB),1):
        #evalue_B[i], evector = np.linalg.eig(NB_U_cart[i])
        #evalue_B[i], evector = LA.eigh(NB_U_cart[i])
        #print(i+1, 'B', evalue_B[i][0]/Angstrom**2, evalue_B[i][1]/Angstrom**2, evalue_B[i][2]/Angstrom**2)
        evalue_B[i] = TE_data_non_proj['thermal_displacements'][0]['displacements'][NB[i]]
        evalue_B[i] = np.sqrt(np.sort(evalue_B[i]))
        eff_disp_B[i] = np.power(evalue_B[i][0]*evalue_B[i][1]*evalue_B[i][2],1./3)
    for i in np.arange(0,len(NO),1):
        #evalue_O[i], evector = np.linalg.eig(NO_U_cart[i])
        #evalue_O[i], evector = LA.eigh(NO_U_cart[i])
        #print(i+1, 'O', evalue_O[i][0]/Angstrom**2, evalue_O[i][1]/Angstrom**2, evalue_O[i][2]/Angstrom**2)
        evalue_O[i] = TE_data_non_proj['thermal_displacements'][0]['displacements'][NO[i]]
        evalue_O[i] = np.sqrt(np.sort(evalue_O[i]))
        eff_disp_O[i] = np.power(evalue_O[i][0]*evalue_O[i][1]*evalue_O[i][2],1./3)
    if persite_amp_True:
        filename = "atom_eff_disp_" + compound_name + ".txt"
        f_atom_eff_disp = open(filename, "w")
        for i in np.arange(0,Natoms,1):
            #evalue_A[i], evector = np.linalg.eig(NA_U_cart[i])
            #evalue_A[i], evector = LA.eigh(NA_U_cart[i])
            #print(i+1, 'A', evalue_A[i][0]/Angstrom**2, evalue_A[i][1]/Angstrom**2, evalue_A[i][2]/Angstrom**2)
            evalue_all_atoms = TE_data_non_proj['thermal_displacements'][0]['displacements'][i]
            evalue_all_atoms = np.sqrt(np.sort(evalue_all_atoms))
            eff_disp_all_atoms = np.power(evalue_all_atoms[0]*evalue_all_atoms[1]*evalue_all_atoms[2],1./3)
            f_atom_eff_disp.write("%d %f\n" %(i+1, eff_disp_all_atoms))

        f_atom_eff_disp.close()

    A_site_TE_mag = np.mean(eff_disp_A)
    B_site_TE_mag = np.mean(eff_disp_B)
    O_site_TE_mag = np.mean(eff_disp_O)

    return A_site_TE_mag, B_site_TE_mag, O_site_TE_mag


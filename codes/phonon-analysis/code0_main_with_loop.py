# KG 7/14/2021

#Todo: AB-first, KPOINTS, 

import sys
import os
import numpy as np
import pandas as pd
import glob
from shutil import copy2
import yaml
import multiprocessing as mp
import re

from code1_POSCARs_info import POSCARs_info_NEB
from code1_POSCARs_info import POSCARs_info_non_NEB
from code2_ev_and_freq_from_phonopy import ev_freq_info_from_phonopy
from code2_ev_and_freq_from_phonopy import phononBC
from code2_ev_and_freq_for_gamma import ev_freq_info_for_gamma
from code3_PBC_correct import PBC_correct
from code4_modal_calc_for_non_gamma import modal_calc
from code4_modal_calc_for_non_gamma import En_hist_n_center
from code4_modal_calc_for_gamma import modal_calc_for_gamma
from code4_modal_calc_for_gamma import En_hist_n_center_for_gamma
import code5_energetics_n_MB_info as Eng_n_MB
import code6_thermal_ellipsoids_for_non_gamma as therm_eps
import code8_thermo_properties as thermo_calc
import code9_bottleneck_area as bottleneck_analysis
import code10_frozen_phonon as frozen_phonon
from code11_OUTCARs_info import OUTCARs_info
import code12_mode_vis as modevis
import code13_seek_path as seekpath
from code14_anharmonic_coupling_constants import do_parallel
from code100_plotting import contributing_phonons_n_PBC_NEB
from code100_plotting import contributing_phonons_n_PBC_non_NEB
import code100_utility_misc as utl_misc

### Parsing the input file
if len(sys.argv) == 1:
    exit("[ERROR] Please, provide the input file...")
inFilename = sys.argv[1]

[Phonopy_path,
DFPT_vasprun_path1,
mode_vis_sortedEn_True,
mode_vis_unsorted_True,
filename_POSCAR_vesta_original,
NEB_True,
NEB_path1,
LEPSILON_calc,
Temperature,
fmin,
MP,
Natoms_in_compound_formula,
FC3_modal_calc,
FC3_file_pathway,
plotting_True,
persite_freq_True,
persite_amp_True] = utl_misc.input_file_parser(inFilename)

#Phonopy_path = "" # If phonopy is not accessible systemwide, give the path to it here
#DFPT_vasprun_path1 = './DFPT_00'   ## folder that includes the DFPT calcs

# Related to mode visualization
#mode_vis_sortedEn_True = False
#mode_vis_sortedEn_True = True
#mode_vis_unsorted_True = False
#mode_vis_unsorted_True = True
#if mode_vis_sortedEn_True or mode_vis_unsorted_True:
    #filename_POSCAR_vesta_original = 'POSCAR.vesta' # Name of the file with all the visualization settings in VESTA for a 39 atom defected cell
    #filename_POSCAR_vesta_original_perfect = 'POSCAR_perfect.vesta' # Name of the file with all the visualization settings in VESTA for a 39 atom defected cell
    #filename_ev_matrix = 'ev.txt'

###
#NEB_True = True
#LEPSILON_calc = False
if NEB_True:
    # Define the images related to your NEB calculation
    image1 = 0 # The image containing the hop origin
    image2 = 1 # The image containing the after the hop image
    image3 = 2 # This will correspond to the TS image number
    image3_S = 2 # for vibrational entropy calculations
    imagef = 4

#Temperature = 300 # K: for thermal ellipsoid and vibrational entropy calculations
#fmin = 0.1 # (THz) minimum frequency included in the PBC, thermal and amplitude calcualtions
# if perovskite:
#Natoms_in_compound_formula = 5 # For ABO3
AB_first = 'A' #### In the perovskite structure: Which atom is declared first in the POSCAR? A-atoms or B-atoms?
#AB_first = 'B'

######
#FC3_modal_calc = True
#FC3_modal_calc = False
#FC3_file_pathway = "./fit3.fcs"
######

###### Start the dataframe
## If the calculations are not based on NEB, there will no value in some of the cells in the dataframe (e.g., NULL)!
data = {
  "Compound": [],
  "Natoms": [],
  "Lattice_volume_A3": [],
  "MB_NCC_eV": [], # from NEB
  "MB_CC_eV": [], # from NEB
  "En_Max_freq_cm": [], # from NEB
  "En_center_freq_cm": [], # from NEB
  "Total_PBC_cm": [],
  "A_PBC_cm": [],
  "B_PBC_cm": [],
  "O_PBC_cm": [],
  "prcnt_neg_freq": [],
  "hopping_O_ave_freq_cm": [], # from NEB
  "hopping_O_eff_disp_A": [], # from NEB
  "hopping_O_proj_eff_disp_A": [], # from NEB
  "hopping_O_max_ellipsoidal_disp_A": [], # from NEB
  "hopping_O_symmetry_index_r1": [], # from NEB
  "hopping_O_symmetry_index_r2": [], # from NEB
  "nearest_A_sites_eff_disp_A": [], # from NEB
  "nearest_B_sites_eff_disp_A": [], # from NEB
  "all_A_sites_eff_disp_ave_A": [],
  "all_B_sites_eff_disp_ave_A": [],
  "all_O_sites_eff_disp_ave_A": [],
  "all_atom_sites_eff_disp ave_A": [],
  "DeltaS_JKmol": [], # from NEB
  "DeltaS_perc_change": [], # from NEB
  "Bottleneck_area_A2": [], # from NEB
  "bn_change_Max_freq_cm": [], # from NEB
  "bn_change_center_freq_cm": [], # from NEB
  "DeltaE_vib_eV": [], # from NEB
  "DeltaE_vib_perc_change": [], # from NEB
  "Dielectric_electronic": [],
  "Dielectric_ionic": [],
  "Rcation_avg_weighted_on_stoichiometry": [],
  "Rdopant": [],
  "RB": [],
  "RA": [],
  "B-O_length_i_2_bonds": [],
  "B-O_length_f_2_bonds": [],
  "A-O_length_i_4_bonds": [],
  "A-O_length_f_4_bonds": [],
  "B_O_length_TS_1_bond": []
}
df = pd.DataFrame(data)

####### Loop over different compounds in the "DFPT_vasprun_path1" folder
f1 = open("whichcompound.txt", "w")
for name in next(os.walk(DFPT_vasprun_path1))[1]:
    compound_name = name
    #if compound_name == 'CaTiO3':
    f1.close()
    f1 = open("whichcompound.txt", "a")
    f1.write(compound_name+'\n')
    f1.close()
    print(compound_name)
    
    # You can ignore the above variables, and just give the folder address directly to the following variable
    DFPT_vasprun_path2 = os.path.join(DFPT_vasprun_path1, compound_name)
    # Our caluclations are based on finished DFPT simulations, so we expect vasprun.xml file to be present in the folder.
    os.system(Phonopy_path + "phonopy --fc " + DFPT_vasprun_path2 + "/vasprun.xml > out.txt") # FORCE_CONSTANTS file should be created in the working directory
    # Hessian should be read here
    # You need to store Hessian in a variable for future use. Otherwise, FORCE CONSTANT file should be read at each q-point loop.

    
    if NEB_True: # Search for NEB folders if we have NEB calculations ready for use
        NEB_path2 = os.path.join(NEB_path1, compound_name)

    # Let's do lattice dynamics using Phonopy: Getting w: freq & v: eigenvectors
    LD_info = ev_freq_info_from_phonopy(DFPT_vasprun_path2, Phonopy_path)
    LD_info.perform_LD_phonopy(MP)

    # Read the structure info (atomic masses, xyz of atoms, etc.)
    if NEB_True:
        struct_info = POSCARs_info_NEB(image1, image2, image3, imagef, NEB_path2) # The passed arguments are the image numbers (folder numbers), TS image number and the address to the NEB folder
        struct_info.read_masses_and_xyzs()

    else:
        struct_info = POSCARs_info_non_NEB(DFPT_vasprun_path2) # The passed argument is the POSCAR address
        struct_info.read_masses_and_xyzs()

    LD_info_for_gamma = ev_freq_info_for_gamma(struct_info.mass) # The passed argument is the atomic mass vector
    LD_info_for_gamma.perform_LD()

    if LEPSILON_calc:
        # Related to reading OUTCAR & dielectric constants
        outcar_info1 = OUTCARs_info(DFPT_vasprun_path2)
        [dielectric_electronic, dielectric_ionic] = outcar_info1.read_dielectric_constants()
    else:
        dielectric_electronic = None
        dielectric_ionic = None

    if NEB_True:
        # calculate diff and apply PBC on it
        dx_direct = struct_info.x1 - struct_info.x2
        diff_correcting = PBC_correct(dx_direct, struct_info.L, struct_info.coeff, struct_info.Natoms) # The passed argument are the diff coords and the L vector & maybe the coeff
        diff_correcting.perform_PBC_correct()

    if NEB_True:
        # Calculate the bottleneck area for perovskite oxide
        [bn_area, NA, NB] = bottleneck_analysis.bottleneck_area_calc(diff_correcting.dx_cartesian,
                                                                     struct_info.x3,
                                                                     struct_info.L,
                                                                     struct_info.coeff,
                                                                     struct_info.Natoms,
                                                                     struct_info.atomNs_passed,
                                                                     struct_info.atomsp_passed)

    if NEB_True:
        # Calculate a few bond lengths for perovskite oxide
        [Rcation_avg_weighted_on_stoichiometry,
         Rdopant,
         RB,
         RA,
         B_O2i, B_O2f, A_O4i, A_O4f, B_O_TS] = bottleneck_analysis.a_few_bond_calc(diff_correcting.dx_cartesian,
                                                                                   struct_info.x1,
                                                                                   struct_info.x3,
                                                                                   struct_info.x4,
                                                                                   struct_info.L,
                                                                                   struct_info.coeff,
                                                                                   struct_info.Natoms,
                                                                                   struct_info.atomNs_passed,
                                                                                   struct_info.atomsp_passed)
    else:
        bn_area = None
        Rcation_avg_weighted_on_stoichiometry = None
        Rdopant = None
        RB = None
        RA = None
        B_O2i = None
        B_O2f = None
        A_O4i = None
        A_O4f = None
        B_O_TS = None

    # Add an underscore before the material_name variable
    compound_name_with_underscore = utl_misc.underscroe_before_numbers(compound_name)

    if NEB_True:
        # Calculate the MB
        MB_NCC = Eng_n_MB.MB_determine(NEB_path2)
        #MB_CC = Eng_n_MB.MB_determine(NEB_path2_forCCMB)
        MB_CC = None
    else:
        MB_NCC = None
        MB_CC = None

    if NEB_True:
        # calculate the modal contribld_inutions
        #modal_info = modal_calc(diff_correcting.dx_cartesian, LD_info.data, struct_info.mass, struct_info.Natoms) # The passed argument are ...
        #modal_info.perform_modal()
        modal_info_for_gamma = modal_calc_for_gamma(diff_correcting.dx_cartesian, LD_info_for_gamma.w, LD_info_for_gamma.v, struct_info.mass, struct_info.Natoms)
        modal_info_for_gamma.perform_modal()

    # Calculate the thermal ellipsoids
    TE_data_non_proj = therm_eps.perform_TE_phonopy(Phonopy_path, fmin, MP)
    if NEB_True: # You need to initialize all the variables to zero and add the contribution calculated for each q-point
        [idx_hopping_ion,
         hopping_ion_TE_mag,
         hopping_ion_TE_projected_mag,
         hopping_ion_TE_max_evalue,
         r1,
         r2] = therm_eps.therm_eps_mag_and_projected_hopping_ion(diff_correcting.dx_cartesian,
                                                                 TE_data_non_proj,
                                                                 struct_info.mass,
                                                                 struct_info.Natoms,
                                                                 Temperature,
                                                                 fmin,
                                                                 Phonopy_path)
        [A_site_TE_mag,
         B_site_TE_mag] = therm_eps.therm_eps_mag_nearest_neighbors_to_hopping_ion(idx_hopping_ion,
                                                                                   struct_info.atomNs_passed,
                                                                                   struct_info.atomsp_passed,
                                                                                   struct_info.x1,
                                                                                   struct_info.L,
                                                                                   struct_info.coeff,
                                                                                   TE_data_non_proj,
                                                                                   struct_info.mass,
                                                                                   struct_info.Natoms,
                                                                                   Temperature,
                                                                                   fmin)
    else:
        hopping_ion_TE_mag = None
        hopping_ion_TE_projected_mag = None
        hopping_ion_TE_max_evalue = None
        A_site_TE_mag = None
        B_site_TE_mag = None
        r1 = None
        r2 = None

    [A_site_TE_mag_all_ave,
     B_site_TE_mag_all_ave,
     O_site_TE_mag_all_ave] = therm_eps.therm_eps_mag_average_all_species(struct_info.atomNs_passed,
                                                                          struct_info.atomsp_passed,
                                                                          struct_info.x1,
                                                                          struct_info.L,
                                                                          struct_info.coeff,
                                                                          TE_data_non_proj,
                                                                          struct_info.mass,
                                                                          struct_info.Natoms,
                                                                          Temperature,
                                                                          fmin,
                                                                          AB_first,
                                                                          compound_name,
                                                                          persite_amp_True)

    all_atom_sites_eff_disp_ave = (1*A_site_TE_mag_all_ave + 1*B_site_TE_mag_all_ave + 3*O_site_TE_mag_all_ave)/5.
    #A_site_TE_mag_all_ave = None
    #B_site_TE_mag_all_ave = None
    #O_site_TE_mag_all_ave = None
    #all_atom_sites_eff_disp_ave = None

    '''if NEB_True:
        # Calculate vibrational entropy
        DFPT_vasprun_path1_image3_S = '../1_TS_vac1_ipv/'
        DFPT_vasprun_path2 = DFPT_vasprun_path1_image3_S + compound_name
        os.system(Phonopy_path + "phonopy --fc " + DFPT_vasprun_path2 + "/vasprun.xml > out.txt") # FORCE_CONSTANTS files should be created in the current folder
        LD_info_image3_S = ev_freq_info(struct_info.mass) # The passed argument is the atomic mass vector
        LD_info_image3_S.perform_LD()
        [diff_vib_S, diff_vib_S_percent_change] = thermo_calc.vibrational_entropy_diff(LD_info.w,
                                                                                       LD_info_image3_S.w,
                                                                                       struct_info.Natoms,
                                                                                       Temperature,
                                                                                       Natoms_in_compound_formula)
    else:
        diff_vib_S = None
        diff_vib_S_percent_change = None

    if NEB_True:
        # Calculate vibrational energy
        DFPT_vasprun_path1_image3_E = '../1_TS_vac1_ipv/'
        DFPT_vasprun_path2 = DFPT_vasprun_path1_image3_E + compound_name
        os.system(Phonopy_path + "phonopy --fc " + DFPT_vasprun_path2 + "/vasprun.xml > out.txt") # FORCE_CONSTANTS files should be created in the current folder
        LD_info_image3_E = ev_freq_info(struct_info.mass) # The passed argument is the atomic mass vector
        LD_info_image3_E.perform_LD()
        [diff_vib_E, diff_vib_E_percent_change] = thermo_calc.vibrational_energy_diff(LD_info.w,
                                                                                      LD_info_image3_E.w,
                                                                                      struct_info.Natoms,
                                                                                      Temperature)
    else:
        diff_vib_E = None
        diff_vib_E_percent_change = None'''

    if NEB_True: # Analysis of the bottleneck area is only possible (for now) when we have a NEB calculation ready
        # Calculate the change in the bottleneck area for perovskite oxide due to modal excitations (frozen phonon analysis)
        # Like En, bn_area_modal can also be a matrix of contributions instead of a vector of contributions
        bn_area_modal = frozen_phonon.modal_bottleneck_area(struct_info.x1,
                                                            struct_info.L,
                                                            struct_info.coeff,
                                                            LD_info_for_gamma.w,
                                                            LD_info_for_gamma.v,
                                                            struct_info.Natoms,
                                                            struct_info.mass,
                                                            NA,
                                                            NB,
                                                            Temperature)
        bn_area_modal_change =  (bn_area_modal - bn_area)/bn_area * 100 # (%)'''

    else:
        bn_area = None
    ### Extract some other numbers
    if NEB_True:
        # Related to modal contributions ############### This section needs to also change in the q-points loop...........
        modalEn_instance = En_hist_n_center_for_gamma(LD_info_for_gamma.w, modal_info_for_gamma.En, 50)
        [w_En, counter_En, w_En_center] = modalEn_instance.hist_for_En_n_center()
        w_En_max = modalEn_instance.freq_max_En()
        #w_En_max *= THz_2_cm
        #w_En_center *= THz_2_cm
    else:
        w_En_max = None
        w_En_center = None
    # Related to Phonon BCs ############### This section needs to also change in the q-points loop...........
    phononBC_instance = phononBC(50, LD_info.data, struct_info.atomNs_passed, struct_info.atomsp_passed, fmin)
    [w_TotalDOS, counter_TotalDOS, TotalPBC] = phononBC_instance.TotalDOS()
    [w_partialDOS, counter_partialDOS, partialPBC] = phononBC_instance.partialDOS('O', AB_first)
    [w_partialDOS_A, counter_partialDOS_A, partialPBC_A] = phononBC_instance.partialDOS('A', AB_first)
    [w_partialDOS_B, counter_partialDOS_B, partialPBC_B] = phononBC_instance.partialDOS('B', AB_first)
    #TotalPBC *= THz_2_cm
    #partialPBC_A *= THz_2_cm
    #partialPBC_B *= THz_2_cm
    #partialPBC *= THz_2_cm

    #partialPBC_A = None
    #partialPBC_B = None
    #partialPBC = None
    if NEB_True:
        # Calculate the average frequency for the hopping O (phonon band center of the hopping O)
        partialPBC_hopping_O = phononBC_instance.ave_freq_hopping_O(idx_hopping_ion)
        #partialPBC_hopping_O *= THz_2_cm
    else:
        partialPBC_hopping_O = None

    if persite_freq_True: # Outputting all the atom-frequencies
        phononBC_instance.all_atom_PBC(compound_name)

    if NEB_True:
        # Related to modal changes in bottleneck (bn) area
        modal_bn_instance = En_hist_n_center_for_gamma(LD_info_for_gamma.w, bn_area_modal_change, 50)
        [w_bn, counter_bn, w_bn_center] = modal_bn_instance.hist_for_En_n_center()
        w_bn_max = modal_bn_instance.freq_max_En()
    else:
        w_bn_max = None
        w_bn_center = None

    # Related to mode visualization
    if mode_vis_sortedEn_True:#################### Put this to be only for q-points = [0, 0, 0] = Gamma-point
        modevis.VESTA_out_POSCAR_sortedEn(diff_correcting.dx_cartesian,
                                          struct_info.x1,
                                          filename_POSCAR_vesta_original,
                                          compound_name,
                                          LD_info_for_gamma.w,
                                          LD_info_for_gamma.v,
                                          modal_info_for_gamma.En)

    if mode_vis_unsorted_True:
        modevis.VESTA_out_POSCAR_unsorted(struct_info.x1,
                                          filename_POSCAR_vesta_original,
                                          compound_name,
                                          LD_info_for_gamma.w,
                                          LD_info_for_gamma.v)

    #for now
    diff_vib_S = None
    diff_vib_S_percent_change = None
    diff_vib_E = None
    diff_vib_E_percent_change = None

    Volume = np.dot(np.cross(struct_info.L[0,:], struct_info.L[1,:]), struct_info.L[2,:])

    # Related to anharmonic coupling constant calculations
    if FC3_modal_calc == True:
        ### Parallel
        Nmodes = struct_info.Natoms * 3
        parallel_class_instance = do_parallel(struct_info.Natoms, struct_info.mass, FC3_file_pathway)
        pool = mp.Pool(mp.cpu_count())
        print(mp.cpu_count())
        results_parallel = pool.map(parallel_class_instance.do_parallel_function, [n for n in np.arange(0, Nmodes, 1)])
        pool.close()

        KmFC3 = np.array(results_parallel)

        filename = "Knml_parallel.txt"
        f_Knml = open(filename, "w")
        for n in np.arange(0, Nmodes, 1):
            for m in np.arange(0, Nmodes, 1):
                for l in np.arange(0, Nmodes, 1):
                    f_Knml.write("%d %d %d %f\n" %(n, m, l, KmFC3[n][m][l]))
        f_Knml.close()

    #### add the info to the dataframe
    df.loc[len(df.index)] = [compound_name_with_underscore, 
                             struct_info.Natoms,
                             Volume,
                             MB_NCC,
                             MB_CC,
                             w_En_max,
                             w_En_center,
                             TotalPBC,
                             partialPBC_A,
                             partialPBC_B,
                             partialPBC,
                             LD_info.prcnt_neg_freq,
                             partialPBC_hopping_O,
                             hopping_ion_TE_mag,
                             hopping_ion_TE_projected_mag,
                             hopping_ion_TE_max_evalue,
                             r1,
                             r2,
                             A_site_TE_mag,
                             B_site_TE_mag,
                             A_site_TE_mag_all_ave,
                             B_site_TE_mag_all_ave,
                             O_site_TE_mag_all_ave,
                             all_atom_sites_eff_disp_ave,
                             diff_vib_S,
                             diff_vib_S_percent_change,
                             bn_area,
                             w_bn_max,
                             w_bn_center,
                             diff_vib_E,
                             diff_vib_E_percent_change,
                             dielectric_electronic,
                             dielectric_ionic,
                             Rcation_avg_weighted_on_stoichiometry,
                             Rdopant,
                             RB,
                             RA,
                             B_O2i,
                             B_O2f,
                             A_O4i,
                             A_O4f,
                             B_O_TS]

    ####Plotting
    if NEB_True and plotting_True: # plot modal contributions and phonon band centers
        plot1 = contributing_phonons_n_PBC_NEB(LD_info_for_gamma.w,
                                               LD_info_for_gamma.v,
                                               modal_info_for_gamma.En,
                                               struct_info.atomsp_passed,
                                               struct_info.atomNs_passed,
                                               compound_name,
                                               bn_area_modal_change,
                                               TotalPBC,
                                               partialPBC,
                                               partialPBC_A,
                                               partialPBC_B,
                                               counter_TotalDOS,
                                               counter_partialDOS,
                                               counter_partialDOS_A,
                                               counter_partialDOS_B,
                                               w_TotalDOS,
                                               w_partialDOS,
                                               w_partialDOS_A,
                                               w_partialDOS_B,)
        plot1.plot_En_PBC(AB_first)
    elif plotting_True: # plot phonon band centers
        '''plot1 = contributing_phonons_n_PBC_non_NEB(LD_info_for_gamma.w,
                                                   LD_info_for_gamma.v,
                                                   struct_info.atomsp_passed,
                                                   struct_info.atomNs_passed,
                                                   compound_name)'''
        plot1 = contributing_phonons_n_PBC_non_NEB(LD_info_for_gamma.w,
                                                  LD_info_for_gamma.v,
                                                  struct_info.atomsp_passed,
                                                  struct_info.atomNs_passed,
                                                  compound_name,
                                                  TotalPBC,
                                                  partialPBC,
                                                  partialPBC_A,
                                                  partialPBC_B,
                                                  counter_TotalDOS,
                                                  counter_partialDOS,
                                                  counter_partialDOS_A,
                                                  counter_partialDOS_B,
                                                  w_TotalDOS,
                                                  w_partialDOS,
                                                  w_partialDOS_A,
                                                  w_partialDOS_B)
        plot1.plot_PBC(AB_first)

    if NEB_True:
        #filename = 'Results_NEB_' + structure_type + '_' + AB_first + '_first_' + str(image1) + '_' + str(image2) + '.csv'
        filename = 'Results_NEB_' + str(image1) + '_' + str(image2) + '.csv'
    else:
        filename = 'Results_nonNEB_' + '.csv'


df.index += 1 
df.to_csv(filename)


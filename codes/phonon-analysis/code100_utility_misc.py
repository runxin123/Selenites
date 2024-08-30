# KG 7/13/2021

import glob
import re

#############################################
def underscroe_before_numbers(val1):
    compound_name_unsplitted = val1
    compound_name_splitted = re.split('(\d+)', compound_name_unsplitted)
    if compound_name_splitted[0] != '':
        compound_name_with_underscore = compound_name_splitted[0]
    counter = 0
    for txt in compound_name_splitted:
        if counter>0 and txt.isalpha() == False and txt != '' and txt != '.' and compound_name_splitted[counter+1] != '.' and compound_name_splitted[counter+1] != '':
            compound_name_with_underscore += txt + '}'
        elif counter>0 and txt.isalpha() == False and txt != '' and txt != '.' and compound_name_splitted[counter+1] == '':
            compound_name_with_underscore += '_' + txt
        elif counter>0 and txt.isalpha() == False and txt != '' and txt != '.' and compound_name_splitted[counter+1] == '.':
            compound_name_with_underscore += '_{' + txt
        elif counter>0 and txt.isalpha() == False and txt != '' and txt == '.':
            compound_name_with_underscore += txt
        elif counter>0 and txt.isalpha() == True and txt != '':
            compound_name_with_underscore += txt
        counter += 1
    return compound_name_with_underscore

############################################# 2/20/2022
def input_file_parser(inFilename):
    with open (inFilename, "r") as infile:
        lines = infile.readlines()

    List_of_needed_input_commands = ['Phonopy_path',
                                     'DFPT_vasprun_path1',
                                     'mode_vis_sortedEn_True',
                                     'mode_vis_unsorted_True',
                                     'filename_POSCAR_vesta_original',
                                     'NEB_True',
                                     'NEB_path1',
                                     'LEPSILON_calc',
                                     'Temperature',
                                     'fmin',
                                     'MP',
                                     'Natoms_in_compound_formula',
                                     'FC3_modal_calc',
                                     'FC3_file_pathway',
                                     'plotting_True',
                                     'persite_freq_True',
                                     'persite_amp_True']

    #Default values
    mode_vis_sortedEn_True = False
    mode_vis_unsorted_True = False
    NEB_True = False
    LEPSILON_calc = False

    for line in lines:
        line = line.replace('\n','')
        line = line.replace('\'','')
        templ = re.split('[=|"|\ ]',line) # splitting the line based on '=', '"', & "\ "
        while("" in templ) :
            templ.remove("")
        if len(templ)>0 and templ[0].split('#')[0]:
            if (templ[0] in List_of_needed_input_commands):
                #Phonopy_path
                if (templ[0] == "Phonopy_path"):
                    if len(templ) == 1:
                        Phonopy_path = ""
                    elif len(templ) > 1:
                        Phonopy_path = templ[1]

                #DFPT_vasprun_path1
                if (templ[0] == "DFPT_vasprun_path1"):
                    if len(templ) == 1:
                        exit("[ERROR] Path to %s not provided... " %(templ[0]))
                    elif len(templ) > 1:
                        DFPT_vasprun_path1 = templ[1]

                #mode_vis_sortedEn_True
                if (templ[0] == "mode_vis_sortedEn_True"):
                    if len(templ) == 1:
                        exit("[ERROR] True or flase shold be mentioned for %s... " %(templ[0]))
                    elif len(templ) > 1:
                        if (templ[1].lower() == 'true'):
                            mode_vis_sortedEn_True = True
                        elif (templ[1].lower() == 'false'):
                            mode_vis_sortedEn_True = False
                        else:
                            exit("[ERROR] Only \"true\" or \"false\" shold be mentioned for %s... " %(templ[0]))

                #mode_vis_unsorted_True
                if (templ[0] == "mode_vis_unsorted_True"):
                    if len(templ) == 1:
                        exit("[ERROR] True or flase shold be mentioned for %s... " %(templ[0]))
                    elif len(templ) > 1:
                        if (templ[1].lower() == 'true'):
                            mode_vis_unsorted_True = True
                        elif (templ[1].lower() == 'false'):
                            mode_vis_unsorted_True = False
                        else:
                            exit("[ERROR] Only \"true\" or \"false\" shold be mentioned for %s... " %(templ[0]))
    
                #filename_POSCAR_vesta_original
                if ((mode_vis_sortedEn_True==True or mode_vis_unsorted_True==True) and templ[0] == "filename_POSCAR_vesta_original"):
                    if len(templ) == 1:
                        exit("[ERROR] Path/name to %s not provided... " %(templ[0]))
                    elif len(templ) > 1:
                        filename_POSCAR_vesta_original = templ[1]
                else:
                    filename_POSCAR_vesta_original = None
    
                #NEB_True
                if (templ[0] == "NEB_True"):
                    if len(templ) == 1:
                        exit("[ERROR] True or flase shold be mentioned for %s... " %(templ[0]))
                    elif len(templ) > 1:
                        if (templ[1].lower() == 'true'):
                            NEB_True = True
                        elif (templ[1].lower() == 'false'):
                            NEB_True = False
                        else:
                            exit("[ERROR] Only \"true\" or \"false\" shold be mentioned for %s... " %(templ[0]))

                #NEB_path1
                if (templ[0] == "NEB_path1"):
                    if len(templ) == 1:
                        exit("[ERROR] Path to %s not provided... " %(templ[0]))
                    elif len(templ) > 1:
                        NEB_path1 = templ[1]
    
                #LEPSILON_calc
                if (templ[0] == "LEPSILON_calc"):
                    if len(templ) == 1:
                        exit("[ERROR] True or flase shold be mentioned for %s... " %(templ[0]))
                    elif len(templ) > 1:
                        if (templ[1].lower() == 'true'):
                            LEPSILON_calc = True
                        elif (templ[1].lower() == 'false'):
                            LEPSILON_calc = False
                        else:
                            exit("[ERROR] Only \"true\" or \"false\" shold be mentioned for %s... " %(templ[0]))
    
                #Temperature
                if (templ[0] == "Temperature"):
                    if len(templ) == 1:
                        exit("[ERROR] value for %s not provided... " %(templ[0]))
                    elif len(templ) > 1:
                        Temperature = float(templ[1])
    
                #fmin
                if (templ[0] == "fmin"):
                    if len(templ) == 1:
                        exit("[ERROR] value for %s not provided... " %(templ[0]))
                    elif len(templ) > 1:
                        fmin = float(templ[1])

                #MP
                if (templ[0] == "MP"):
                    if len(templ) == 1:
                        exit("[ERROR] value for %s not provided... " %(templ[0]))
                    elif len(templ) > 1:
                        MP = [templ[1] + ' ' + templ[1] + ' ' + templ[1]]

                #Natoms_in_compound_formula
                if (templ[0] == "Natoms_in_compound_formula"):
                    if len(templ) == 1:
                        exit("[ERROR] value for %s not provided... " %(templ[0]))
                    elif len(templ) > 1:
                        Natoms_in_compound_formula = int(templ[1])
    
                #FC3_modal_calc
                if (templ[0] == "FC3_modal_calc"):
                    if len(templ) == 1:
                        exit("[ERROR] True or flase shold be mentioned for %s... " %(templ[0]))
                    elif len(templ) > 1:
                        if (templ[1].lower() == 'true'):
                            FC3_modal_calc = True
                        elif (templ[1].lower() == 'false'):
                            FC3_modal_calc = False
                        else:
                            exit("[ERROR] Only \"true\" or \"false\" shold be mentioned for %s... " %(templ[0]))
    
                #FC3_file_pathway
                if (templ[0] == "FC3_file_pathway"):
                    if len(templ) == 1:
                        exit("[ERROR] Path to %s not provided... " %(templ[0]))
                    elif len(templ) > 1:
                        FC3_file_pathway = templ[1]

                #plotting_True
                if (templ[0] == "plotting_True"):
                    if len(templ) == 1:
                        exit("[ERROR] True or flase shold be mentioned for %s... " %(templ[0]))
                    elif len(templ) > 1:
                        if (templ[1].lower() == 'true'):
                            plotting_True = True
                        elif (templ[1].lower() == 'false'):
                            plotting_True = False
                        else:
                            exit("[ERROR] Only \"true\" or \"false\" shold be mentioned for %s... " %(templ[0]))

                #persite_freq_True
                if (templ[0] == "persite_freq_True"):
                    if len(templ) == 1:
                        exit("[ERROR] True or flase shold be mentioned for %s... " %(templ[0]))
                    elif len(templ) > 1:
                        if (templ[1].lower() == 'true'):
                            persite_freq_True = True
                        elif (templ[1].lower() == 'false'):
                            persite_freq_True = False
                        else:
                            exit("[ERROR] Only \"true\" or \"false\" shold be mentioned for %s... " %(templ[0]))

                #persite_amp_True
                if (templ[0] == "persite_amp_True"):
                    if len(templ) == 1:
                        exit("[ERROR] True or flase shold be mentioned for %s... " %(templ[0]))
                    elif len(templ) > 1:
                        if (templ[1].lower() == 'true'):
                            persite_amp_True = True
                        elif (templ[1].lower() == 'false'):
                            persite_amp_True = False
                        else:
                            exit("[ERROR] Only \"true\" or \"false\" shold be mentioned for %s... " %(templ[0]))
    
            else:
                exit("[ERROR] Command \"%s\" does not exist." %(templ[0]))        

    return Phonopy_path, DFPT_vasprun_path1,\
           mode_vis_sortedEn_True,\
           mode_vis_unsorted_True,\
           filename_POSCAR_vesta_original,\
           NEB_True,\
           NEB_path1,\
           LEPSILON_calc,\
           Temperature,\
           fmin,\
           MP[0],\
           Natoms_in_compound_formula,\
           FC3_modal_calc,\
           FC3_file_pathway,\
           plotting_True,\
           persite_freq_True,\
           persite_amp_True


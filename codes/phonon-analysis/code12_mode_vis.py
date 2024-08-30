# KG 10/19/2021

import numpy as np
import shutil
from pathlib import Path

def VESTA_out_POSCAR_sortedEn(dx, x_coords, filename_POSCAR_original,filename1, compound_name, w, v, En):
    # Some initital parameters
    scale = 8
    radius = 0.25
    RGB = "255 0 0"

    Natoms = np.size(x_coords,0)

    # find the index of the hopping ion
    maxdx_idx = 0
    maxdx = 0
    for i in np.arange(0,Natoms,1):
        temp = np.sqrt(np.sum(np.square(dx[i])))
        if temp > maxdx:
            maxdx_idx = i
            maxdx = temp
    idx_hopping_ion = maxdx_idx
    #print(idx_hopping_ion)

    # Read the original POSCAR_vesta file
    with open(filename_POSCAR_original, "r") as f:
        lines = f.readlines()

    # Create the mode_vis directory
    
    filepath1 = "./" + filename1
    #try:
    #    shutil.rmtree(filepath1)
    #except OSError as e:
    #    pass
    #Path(filepath1).mkdir(parents=True, exist_ok=True)

    # Create the compound_name directory
    filename2 = compound_name + '/'
    filepath2 = filepath1 + filename2
    try:
        shutil.rmtree(filepath2)
    except OSError as e:
        pass
    Path(filepath2).mkdir(parents=True, exist_ok=True)

    # Sort En and get the sorted arguments
    args = np.argsort(En) # args sorted from min to max
    args = args[::-1] # args sorted from max to min


    for n in np.arange(0, len(args), 1): # This loops over all the modes, but when they are sorted based on their En
        # Open the Vesta file for visualizing the mode n ()n index for the degree of contribution
        ev = v[:,args[n]]
        freq_str = ('%.2f' %(w[args[n]]))
        filename_POSCARvesta = 'POSCAR_' + compound_name + '_c' + str(n+1) + '_n' + str(args[n]+1) + '_f' + freq_str + '.vesta'
        f_poscar_temp = open(filepath2 + filename_POSCARvesta, "w")
        
        countline = 0
        coord_done = 0
        while countline < len(lines):
            line_splitted = lines[countline].split()
            if line_splitted and coord_done == 0 and line_splitted[0] == 'STRUC': # Insert the coordinates in this section
                f_poscar_temp.write(lines[countline])
                countline += 1
                for i in range(0, Natoms):
                    line_splitted = lines[countline].split()
                    line_splitted[4] = ('%.6f' %(x_coords[i][0]))
                    line_splitted[5] = ('%.6f' %(x_coords[i][1]))
                    line_splitted[6] = ('%.6f' %(x_coords[i][2]))
                    lines[countline] = ''
                    for j in np.arange(0, len(line_splitted), 1):
                        lines[countline] += line_splitted[j] + ' '
                    lines[countline] += '\n'
                    f_poscar_temp.write(lines[countline])
                    countline += 1
                    f_poscar_temp.write(lines[countline])
                    countline += 1
                coord_done = 1


            elif line_splitted and line_splitted[0] == 'SITET': # Change the color for the hopping ion to yellow!
                f_poscar_temp.write(lines[countline])
                countline += 1
                for i in range(0, Natoms):
                    if i == idx_hopping_ion:
                        line_splitted = lines[countline].split()
                        line_splitted[3] = ('255') # R
                        line_splitted[4] = ('255') # G
                        line_splitted[5] = ('0') # B
                        lines[countline] = ''
                        for j in np.arange(0, len(line_splitted), 1):
                            lines[countline] += line_splitted[j] + ' '
                        lines[countline] += '\n'
                        f_poscar_temp.write(lines[countline])
                        countline += 1
                    else:
                        f_poscar_temp.write(lines[countline])
                        countline += 1
            elif line_splitted and line_splitted[0] == 'VECTR': # Insert the eigenvectors in this section
                f_poscar_temp.write(lines[countline])
                countline += 1
                for i in range(0, Natoms):
                    mystring = str(i+1) + ' '
                    mystring += ('%.6f' %(scale*ev[3*i+0])) + ' '
                    mystring += ('%.6f' %(scale*ev[3*i+1])) + ' '
                    mystring += ('%.6f' %(scale*ev[3*i+2])) + ' '
                    mystring += '0\n'
                    f_poscar_temp.write(mystring)
                    mystring = str(i+1) + ' 0 0 0 0' + '\n' + '0 0 0 0 0' + '\n'
                    f_poscar_temp.write(mystring)
                f_poscar_temp.write(lines[countline]) # add the last "0 0 0 0 0" at the end of "VECTR" section
                countline += 1
            elif line_splitted and line_splitted[0] == 'VECTT': # Insert radius & RGB data in this section
                f_poscar_temp.write(lines[countline])
                countline += 1
                for i in range(0, Natoms):
                    mystring = str(i+1) + ' ' + str(radius) + ' ' + RGB + ' 0\n'
                    f_poscar_temp.write(mystring)
                f_poscar_temp.write(lines[countline]) # add the last "0 0 0 0 0" at the end of "VECTT" section
                countline += 1
            else:
                f_poscar_temp.write(lines[countline])
                countline += 1
        f_poscar_temp.close()

def VESTA_out_POSCAR_unsorted(x_coords, filename_POSCAR_original, filename1,compound_name, w, v):
    # Some initital parameters
    scale = 8
    radius = 0.25
    RGB = "255 0 0"

    Natoms = np.size(x_coords,0)
    Nmodes = Natoms * 3


    # Read the original POSCAR_vesta file
    with open(filename_POSCAR_original, "r") as f:
        lines = f.readlines()

    # Create the mode_vis directory
    # filename1 = 'mode_vis_sorted_based_on_mode_id/'
    filepath1 = "./" + filename1
    #try:
    #    shutil.rmtree(filepath1)
    #except OSError as e:
    #    pass
    #Path(filepath1).mkdir(parents=True, exist_ok=True)

    # Create the compound_name directory
    filename2 = compound_name + '/'
    filepath2 = filepath1 + filename2
    try:
        shutil.rmtree(filepath2)
    except OSError as e:
        pass
    Path(filepath2).mkdir(parents=True, exist_ok=True)

    for n in np.arange(0, Nmodes, 1):
        ev = v[:,n]
        freq_str = ('%.2f' %(w[n]))
        filename_POSCARvesta = 'POSCAR_' + compound_name + '_n' + str(n+1) + '_f' + freq_str + '.vesta'
        f_poscar_temp = open(filepath2 + filename_POSCARvesta, "w")
        
        countline = 0
        coord_done = 0
        while countline < len(lines):
            line_splitted = lines[countline].split()
            if line_splitted and coord_done == 0 and line_splitted[0] == 'STRUC': # Insert the coordinates in this section
                f_poscar_temp.write(lines[countline])
                countline += 1
                for i in range(0, Natoms):
                    line_splitted = lines[countline].split()
                    line_splitted[4] = ('%.6f' %(x_coords[i][0]))
                    line_splitted[5] = ('%.6f' %(x_coords[i][1]))
                    line_splitted[6] = ('%.6f' %(x_coords[i][2]))
                    lines[countline] = ''
                    for j in np.arange(0, len(line_splitted), 1):
                        lines[countline] += line_splitted[j] + ' '
                    lines[countline] += '\n'
                    f_poscar_temp.write(lines[countline])
                    countline += 1
                    f_poscar_temp.write(lines[countline])
                    countline += 1
                coord_done = 1
            elif line_splitted and line_splitted[0] == 'SITET':
                f_poscar_temp.write(lines[countline])
                countline += 1
                for i in range(0, Natoms):
                    f_poscar_temp.write(lines[countline])
                    countline += 1
            elif line_splitted and line_splitted[0] == 'VECTR': # Insert the eigenvectors in this section
                f_poscar_temp.write(lines[countline])
                countline += 1
                for i in range(0, Natoms):
                    mystring = str(i+1) + ' '
                    mystring += ('%.6f' %(scale*ev[3*i+0])) + ' '
                    mystring += ('%.6f' %(scale*ev[3*i+1])) + ' '
                    mystring += ('%.6f' %(scale*ev[3*i+2])) + ' '
                    mystring += '0\n'
                    f_poscar_temp.write(mystring)
                    mystring = str(i+1) + ' 0 0 0 0' + '\n' + '0 0 0 0 0' + '\n'
                    f_poscar_temp.write(mystring)
                f_poscar_temp.write(lines[countline]) # add the last "0 0 0 0 0" at the end of "VECTR" section
                countline += 1
            elif line_splitted and line_splitted[0] == 'VECTT': # Insert radius & RGB data in this section
                f_poscar_temp.write(lines[countline])
                countline += 1
                for i in range(0, Natoms):
                    mystring = str(i+1) + ' ' + str(radius) + ' ' + RGB + ' 0\n'
                    f_poscar_temp.write(mystring)
                f_poscar_temp.write(lines[countline]) # add the last "0 0 0 0 0" at the end of "VECTT" section
                countline += 1
            else:
                f_poscar_temp.write(lines[countline])
                countline += 1
        f_poscar_temp.close()


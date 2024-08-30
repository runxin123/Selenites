# KG 7/13/2021

import os
import numpy as np

def MB_determine(val1):
    NEB_path = val1
    list_of_folders = next(os.walk(NEB_path))[1]
    Num_images = len(list_of_folders)-1 # Number of images (folders)
    E = np.zeros(Num_images)
    for i in np.arange(0,Num_images,1):
        Efilename = os.path.join(NEB_path, '0' + str(i))
        with open(os.path.join(Efilename, "OSZICAR"), "r") as f:
            lines = f.readlines()
            templ = lines[len(lines)-1].split()
            E[i] = float(templ[4])
    max_deltaE = 0
    for i in np.arange(0,Num_images,1):
        for j in np.arange(i,Num_images,1):
            temp = abs(E[i]-E[j])
            if temp > max_deltaE:
                max_deltaE = temp
    return max_deltaE

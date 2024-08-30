# KG 10/13/2021

import numpy as np

class OUTCARs_info:
    def __init__(self, val1):
        self.path_OUTCAR = val1

    def read_dielectric_constants(self):
        macroscopic_static_dielectric_tensor = np.zeros([3, 3])
        macroscopic_dielectric_tensor_ionic_contribution = np.zeros([3, 3]) 
        filename = self.path_OUTCAR + "/OUTCAR"
        with open(filename, "r") as f:
            counter = 0
            found_light1 = 0
            lineafter1 = 0
            found_light2 = 0
            lineafter2 = 0
            break1 = 0
            break2 = 0
            first_tensor = 0 # This is to avoid reading the "MACROSCOPIC STATIC DIELECTRIC TENSOR" two times, since we have two instances of it in the OUTCAR.
            for line in f:
                if len(line.split()) > 0:
                    '''# Get the "HEAD OF MICROSCOPIC STATIC DIELECTRIC TENSOR"
                    if line.split()[0] == 'HEAD':
                        found_light1 = 1
                    if found_light1 == 1 and lineafter1 > 1 and lineafter1 < 5:
                        macroscopic_dielectric_tensor[lineafter1-2][0] = float(line.split()[0])
                        macroscopic_dielectric_tensor[lineafter1-2][1] = float(line.split()[1])
                        macroscopic_dielectric_tensor[lineafter1-2][2] = float(line.split()[2])
                        lineafter1 += 1
                        if lineafter1 == 5:
                            found_light1 = 0
                            break1 = 1
                            if break1 == 1 and break2 == 1:
                                break
                    elif found_light1 == 1:  # Just to go past the '------' line
                        lineafter1 += 1'''

                    # Get the "MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)"
                    if line.split()[0] == 'MACROSCOPIC' and \
                       line.split()[1] == 'STATIC' and \
                       line.split()[2] == 'DIELECTRIC' and \
                       line.split()[3] == 'TENSOR' and \
                       line.split()[4] == '(including' and \
                       first_tensor == 0:
                        #print(line.split())
                        first_tensor = 1
                        found_light1 = 1
                    if found_light1 == 1 and lineafter1 > 1 and lineafter1 < 5:
                        macroscopic_static_dielectric_tensor[lineafter1-2][0] = float(line.split()[0])
                        macroscopic_static_dielectric_tensor[lineafter1-2][1] = float(line.split()[1])
                        macroscopic_static_dielectric_tensor[lineafter1-2][2] = float(line.split()[2])
                        lineafter1 += 1
                        if lineafter1 == 5:
                            found_light1 = 0
                            break1 = 1
                            if break1 == 1 and break2 == 1:
                                break
                    elif found_light1 == 1:  # Just to go past the '------' line
                        lineafter1 += 1

                    # Get the " MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION"
                    if line.split()[0] == 'MACROSCOPIC' and line.split()[4] == 'IONIC':
                        found_light2 = 1
                    if found_light2 == 1 and lineafter2 > 1 and lineafter2 < 5:
                        macroscopic_dielectric_tensor_ionic_contribution[lineafter2-2][0] = float(line.split()[0])
                        macroscopic_dielectric_tensor_ionic_contribution[lineafter2-2][1] = float(line.split()[1])
                        macroscopic_dielectric_tensor_ionic_contribution[lineafter2-2][2] = float(line.split()[2])
                        lineafter2 += 1
                        if lineafter2 == 5:
                            found_light2 = 0
                            break2 = 1
                            if break1 == 1 and break2 == 1:
                                break
                    elif found_light2 == 1:  # Just to go past the '------' line
                        lineafter2 += 1

                counter += 1
        dielectric_electronic = (macroscopic_static_dielectric_tensor[0][0] + \
                                 macroscopic_static_dielectric_tensor[1][1] + \
                                 macroscopic_static_dielectric_tensor[2][2])/3
        dielectric_ionic = (macroscopic_dielectric_tensor_ionic_contribution[0][0] + \
                            macroscopic_dielectric_tensor_ionic_contribution[1][1] + \
                            macroscopic_dielectric_tensor_ionic_contribution[2][2])/3

        #print(macroscopic_static_dielectric_tensor)

        return dielectric_electronic, dielectric_ionic


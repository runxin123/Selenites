import numpy as np
import multiprocessing as mp
import re

class do_parallel:
    def __init__(self, val1, val2, val3):
        self.Natoms = val1
        self.mass=val2
        Natoms = self.Natoms
        Nmodes = Natoms * 3
        filename = val3

        self.FC3 = np.zeros([Natoms, Natoms, Natoms, 3, 3, 3])
        self.ev = np.zeros([Nmodes, Nmodes])
        #################### Reading FC3 values
        with open(filename, "r") as f:
            lines = f.readlines()

        counter = 0
        while True:
            templ = lines[counter].split()
            if len(templ) > 0:
                if templ[0] == '*FC3':
                    break
            counter += 1

        NFC3 = 0
        while True:
            counter += 1
            templ = lines[counter].split()
            if len(templ) > 0:
                FC = float(templ[2])
                atomi = int(re.split('(\d+)',templ[4])[1])-1
                atomj = int(re.split('(\d+)',templ[5])[1])-1
                atomk = int(re.split('(\d+)',templ[6])[1])-1
                dira = self.find_dir(re.split('(\d+)',templ[4])[2])
                dirb = self.find_dir(re.split('(\d+)',templ[5])[2])
                dirg = self.find_dir(re.split('(\d+)',templ[6])[2])
                self.FC3[atomi][atomj][atomk][dira][dirb][dirg] = FC
                NFC3 += 1
            else:
                break

        #################### Reading eigen-vector values
        filename = "./ev.txt"
        with open(filename, "r") as f:
            lines = f.readlines()

        if len(lines) != Nmodes:
            print("Something is wrong!!!!!!!")

        for i in np.arange(0, Nmodes, 1):
            templ = lines[i].split()
            for j in np.arange(0, Nmodes, 1):
                self.ev[i][j] = float(templ[j])

    def find_dir(self, dir_letter):
        if dir_letter == 'x':
            return 0
        elif dir_letter == 'y':
            return 1
        elif dir_letter == 'z':
            return 2
    
    def do_serial_function(self):
        Natoms = self.Natoms # This should be replaced by copy()
        Nmodes = Natoms * 3
        mass = self.mass # This should be replaced by copy()
        FC3 = self.FC3
        ev = self.ev
        KmFC3 = np.zeros([Nmodes, Nmodes, Nmodes]) #K modes FC3

        #################### Calculating the mode coupling constants K_nml
        for n in np.arange(0, Nmodes, 1):
            print(n)
            for m in np.arange(0, Nmodes, 1):
                for l in np.arange(0, Nmodes, 1):
                    sum1 = 0
                    for i in np.arange(0, Natoms, 1):
                        for j in np.arange(0, Natoms, 1):
                            for k in np.arange(0, Natoms, 1):
                                for a in np.arange(0, 3, 1): # alpha
                                    for b in np.arange(0, 3, 1): # beta
                                        for g in np.arange(0, 3, 1): # gamma
                                            #print(i)
                                            if FC3[i][j][k][a][b][g] != 0:
                                                sum1 += FC3[i][j][k][a][b][g]*ev[i*3+a][n]*ev[j*3+b][m]*ev[k*3+g][l]/(np.sqrt(mass[i]*mass[j]*mass[k]))
                    KmFC3[n][m][l] = sum1

        filename = "Knml.txt"
        f_Knml = open(filename, "w")
        for n in np.arange(0, Nmodes, 1):
            for m in np.arange(0, Nmodes, 1):
                for l in np.arange(0, Nmodes, 1):
                    f_Knml.write("%d %d %d %f\n" %(n, m, l, KmFC3[n][m][l]))
        f_Knml.close()
        print(KmFC3.shape)

    def do_parallel_function(self, n):
        print(n) # Counter for the progress of the parallel processing
        f2 = open("whichcompound_progress.txt", "a")
        f2.write(str(n)+'\n')
        f2.close()
        Natoms = self.Natoms
        Nmodes = Natoms * 3
        mass = self.mass # This should be replaced by copy()
        FC3 = self.FC3
        ev = self.ev

        #################### Calculating the mode coupling constants K_nml
        results_temp = np.zeros([Nmodes, Nmodes])
        for m in np.arange(0, Nmodes, 1):
            for l in np.arange(0, Nmodes, 1):
                sum1 = 0
                for i in np.arange(0, Natoms, 1):
                    for j in np.arange(0, Natoms, 1):
                        for k in np.arange(0, Natoms, 1):
                            for a in np.arange(0, 3, 1): # alpha
                                for b in np.arange(0, 3, 1): # beta
                                    for g in np.arange(0, 3, 1): # gamma
                                        if FC3[i][j][k][a][b][g] != 0:
                                            sum1 += FC3[i][j][k][a][b][g]*ev[i*3+a][n]*ev[j*3+b][m]*ev[k*3+g][l]/(np.sqrt(mass[i]*mass[j]*mass[k]))
                results_temp[m][l] = sum1

        return results_temp


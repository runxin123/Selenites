# KG 11/8/2021

import numpy as np
from shutil import copy2
import yaml
import os

class ev_freq_info_from_phonopy:
    def __init__(self, val1, val2):
        self.DFPT_vasprun_path2 = val1
        self.Phonopy_path = val2

    def perform_LD_phonopy(self, MP,): # Function to calculate frequencies and eigenvectors at all q-points
        # copy2(self.DFPT_vasprun_path2 + '/POSCAR',)
        # os.chdir(filename1)
        # copy2('../FORCE_CONSTANTS','.')

        # with open("./mesh_0.conf", "r") as f:   ####
        #     lines = f.readlines()
        # for i in range(0, len(lines)):
        #     templ = lines[i].split()
        #     if templ[0].lower() == 'mp':
        #         lines[i] = 'MP = ' + MP + '\n'
        # filename = "/lustre/home/acct-umjzhh/umjzhh-5/oyrx/MN-rules/mesh.conf"
        # f_mesh = open(filename, "w")
        # for i in range(0, len(lines)):
        #     f_mesh.write("%s\n" %(lines[i][:-1]))
        # f_mesh.close()

        # os.system(self.Phonopy_path + "phonopy --readfc --dos mesh.conf > out.txt")
        # os.chdir('..')

        #print('I am here!!!!')
        with open(f"{self.Phonopy_path}mesh.yaml") as f:
            self.data = yaml.load(f, Loader=yaml.FullLoader)
        #print('I am here2!!!!')
        # Let's calculate the percentage of negative frequencies
        count_total = 0
        count_neg_freq = 0
        N_qpoints = len(self.data['phonon'])
        N_modes_each_qpoints = len(self.data['phonon'][0]['band'])
        for nq in np.arange(0, N_qpoints, 1):
            w = np.zeros([N_modes_each_qpoints]) # Initiating the freq array
            for nm in np.arange(0, N_modes_each_qpoints, 1):
                w[nm] = self.data['phonon'][nq]['band'][nm]['frequency']
                count_total += 1
                if w[nm] < 0:
                    count_neg_freq += 1
        self.prcnt_neg_freq = count_neg_freq/count_total*100

class phononBC:
    def __init__(self,val1, val2, val3, val4, val5):
        self.Npts = val1
        self.data = val2
        self.atomNs = val3
        self.atomsp = val4
        self.fmin = val5

    def TotalDOS(self): # Function to calculate Total phonon BCs
        N_qpoints = len(self.data['phonon'])
        N_modes_each_qpoints = len(self.data['phonon'][0]['band']) # Think of this as Nband = Nmodes
        #if (N_modes_each_qpoints != Nmodes):
        #    print('Waaaaaaarnning!!!!!!!!!!!!')
        Natoms = int(N_modes_each_qpoints/3)
        counter = np.zeros(self.Npts) # Histogram bins
        # Let's put all the frequencies in a big vector
        freq = np.zeros([N_qpoints*N_modes_each_qpoints]) # Initiating the freq array
        countfreq = 0
        for nq in np.arange(0, N_qpoints, 1):
            for nm in np.arange(0, N_modes_each_qpoints, 1):
                freq[countfreq] = self.data['phonon'][nq]['band'][nm]['frequency']
                countfreq += 1
        min_freq = min(freq)
        max_freq = max(freq)
        domega = (max_freq - min_freq)/self.Npts
        w = np.arange(1,self.Npts+1,1)*domega + min_freq + domega*.5

        # Calculate the total DOS
        for nq in np.arange(0, N_qpoints, 1):
            # Let's prepare the freq vector and ev matrix for each q-point first 
            freq = np.zeros([N_modes_each_qpoints]) # Initiating the freq array
            v = np.zeros([N_modes_each_qpoints, N_modes_each_qpoints], dtype=complex) # Initiating the ev array
            for nm in np.arange(0, N_modes_each_qpoints, 1): # This loop is not needed. It will be done inside each q-point calculation.
                freq[nm] = self.data['phonon'][nq]['band'][nm]['frequency']
            for i in np.arange(0,np.size(freq),1):
                if (np.floor((freq[i]-min_freq)/domega)<self.Npts):
                    indx = np.floor((freq[i]-min_freq)/domega)
                else:
                    indx = np.floor((freq[i]-min_freq)/domega)-1
                counter[int(indx)] += 1
        # Calculate the total PBC
        sum1 = 0
        sum2 = 0
        for i in np.arange(0,np.size(counter),1):
            if (w[i] > self.fmin):
                sum1 += w[i] * counter[i]
                sum2 += counter[i]
        totalPBC = sum1/sum2
        return w, counter, totalPBC

    def partialDOS(self, val1, val2): # Function to calculate partial phonon BCs
        spcs_name = val1
        AB_first = val2
        myatomNS = self.atomNs.copy()
        # Calculate the correct range for performing the summation
        if (AB_first == 'A'):
            if len(self.atomsp) == 3:
                if spcs_name == 'A':
                    spcs_name = self.atomsp[0]
                elif spcs_name == 'B':
                    spcs_name = self.atomsp[1]
                counter = 0
                for spname in self.atomsp:
                    if spname == spcs_name:
                        break
                    counter += 1
            elif len(self.atomsp) == 4:
                if spcs_name == 'A':
                    spcs_name = self.atomsp[1] # This is based on the preparation of all the POSCARs
                    myatomNS[1] = myatomNS[0]+myatomNS[1]
                    myatomNS[0] = 0
                elif spcs_name == 'B':
                    spcs_name = self.atomsp[2]
                counter = 0
                for spname in self.atomsp:
                    if spname == spcs_name:
                        break
                    counter += 1

            elif len(self.atomsp) == 5:
                if spcs_name == 'A':
                    spcs_name = self.atomsp[1] # This is based on the preparation of all the POSCARs
                    myatomNS[1] = myatomNS[0]+myatomNS[1]
                    myatomNS[0] = 0
                elif spcs_name == 'B':
                    spcs_name = self.atomsp[2]
                counter = 0
                for spname in self.atomsp:
                    if spname == spcs_name:
                        break
                    counter += 1

        elif (AB_first == 'B'):
            if len(self.atomsp) == 3:
                if spcs_name == 'A':
                    spcs_name = self.atomsp[1]
                elif spcs_name == 'B':
                    spcs_name = self.atomsp[0]
                counter = 0
                for spname in self.atomsp:
                    if spname == spcs_name:
                        break
                    counter += 1
            elif len(self.atomsp) == 4:
                if spcs_name == 'A':
                    spcs_name = self.atomsp[2] # This is based on the preparation of all the POSCARs
                    myatomNS[2] = self.atomNs[1]+self.atomNs[2]
                    myatomNS[1] = self.atomNs[0]
                    myatomNS[0] = 0
                elif spcs_name == 'B':
                    spcs_name = self.atomsp[0]
                counter = 0
                for spname in self.atomsp:
                    if spname == spcs_name:
                        break
                    counter += 1
        begidx = 0
        endidx = 0
        '''for i in np.arange(0,counter+1,1):
            begidx = endidx
            endidx += myatomNS[i]'''
        begidx = self.atomNs[0]-1
        endidx = begidx + self.atomNs[1]
        #print(self.atomNs)
        #print('%d %d' %(begidx, endidx))

        N_qpoints = len(self.data['phonon'])
        N_modes_each_qpoints = len(self.data['phonon'][0]['band']) # Think of this as Nband = Nmodes
        #if (N_modes_each_qpoints != Nmodes):
        #    print('Waaaaaaarnning!!!!!!!!!!!!')
        Natoms = int(N_modes_each_qpoints/3)
        #print(Natoms)
        counter = np.zeros(self.Npts) # Histogram bins
        # Let's put all the frequencies in a big vector
        freq = np.zeros([N_qpoints*N_modes_each_qpoints]) # Initiating the freq array
        countfreq = 0
        for nq in np.arange(0, N_qpoints, 1):
            for nm in np.arange(0, N_modes_each_qpoints, 1):
                freq[countfreq] = self.data['phonon'][nq]['band'][nm]['frequency']
                countfreq += 1
        min_freq = min(freq)
        max_freq = max(freq)
        domega = (max_freq - min_freq)/self.Npts
        w = np.arange(1,self.Npts+1,1)*domega + min_freq + domega*.5

        # Calculate the total DOS
        for nq in np.arange(0, N_qpoints, 1):
            # Let's prepare the freq vector and ev matrix for each q-point first 
            freq = np.zeros([N_modes_each_qpoints]) # Initiating the freq array
            v = np.zeros([N_modes_each_qpoints, N_modes_each_qpoints], dtype=complex) # Initiating the ev array
            for nm in np.arange(0, N_modes_each_qpoints, 1): # This loop is not needed. It will be done inside each q-point calculation.
                freq[nm] = self.data['phonon'][nq]['band'][nm]['frequency']
                for na in np.arange(0, Natoms, 1):
                    v.real[na*3+0, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][0][0]
                    v.real[na*3+1, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][1][0]
                    v.real[na*3+2, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][2][0]
                    v.imag[na*3+0, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][0][1]
                    v.imag[na*3+1, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][1][1]
                    v.imag[na*3+2, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][2][1]

            # Calculate the partial DOS - it is the summation of the effect coming from different q-points
            for i in np.arange(0,np.size(freq),1):
                if (np.floor((freq[i]-min_freq)/domega)<self.Npts):
                    indx = np.floor((freq[i]-min_freq)/domega)
                else:
                    indx = np.floor((freq[i]-min_freq)/domega)-1
                sum1 = 0
                for j in np.arange(begidx,endidx,1):
                    for k in np.arange(0,3):
                        #sum1 += np.power(v[3*j+k][i],2)
                        sum1 += np.power(np.absolute(v[3*j+k][i]),2)
                counter[int(indx)] += sum1
        # Calculate the partial PBC
        sum1 = 0
        sum2 = 0
        for i in np.arange(0,np.size(counter),1):
            if (w[i] > self.fmin):
                sum1 += w[i] * counter[i]
                sum2 += counter[i]
        partialPBC = sum1/sum2
        return w, counter, partialPBC

    def ave_freq_hopping_O(self, val1): # Function to calculate the average frequency of the hopping ion
        idx_hopping_ion = val1
        print("Index hopping atom: %d" %(idx_hopping_ion+1))
        N_qpoints = len(self.data['phonon'])
        N_modes_each_qpoints = len(self.data['phonon'][0]['band']) # Think of this as Nband = Nmodes
        #if (N_modes_each_qpoints != Nmodes):
        #    print('Waaaaaaarnning!!!!!!!!!!!!')
        Natoms = int(N_modes_each_qpoints/3)
        counter = np.zeros(self.Npts) # Histogram bins
        # Let's put all the frequencies in a big vector
        freq = np.zeros([N_qpoints*N_modes_each_qpoints]) # Initiating the freq array
        countfreq = 0
        for nq in np.arange(0, N_qpoints, 1):
            for nm in np.arange(0, N_modes_each_qpoints, 1):
                freq[countfreq] = self.data['phonon'][nq]['band'][nm]['frequency']
                countfreq += 1
        min_freq = min(freq)
        max_freq = max(freq)
        domega = (max_freq - min_freq)/self.Npts
        w = np.arange(1,self.Npts+1,1)*domega + min_freq + domega*.5

        for nq in np.arange(0, N_qpoints, 1):
            # Let's prepare the freq vector and ev matrix for each q-point first 
            freq = np.zeros([N_modes_each_qpoints]) # Initiating the freq array
            v = np.zeros([N_modes_each_qpoints, N_modes_each_qpoints], dtype=complex) # Initiating the ev array
            for nm in np.arange(0, N_modes_each_qpoints, 1): # This loop is not needed. It will be done inside each q-point calculation.
                freq[nm] = self.data['phonon'][nq]['band'][nm]['frequency']
                for na in np.arange(0, Natoms, 1):
                    v.real[na*3+0, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][0][0]
                    v.real[na*3+1, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][1][0]
                    v.real[na*3+2, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][2][0]
                    v.imag[na*3+0, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][0][1]
                    v.imag[na*3+1, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][1][1]
                    v.imag[na*3+2, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][2][1]

            for i in np.arange(0,np.size(freq),1):
                if (np.floor((freq[i]-min_freq)/domega)<self.Npts):
                    indx = np.floor((freq[i]-min_freq)/domega)
                else:
                    indx = np.floor((freq[i]-min_freq)/domega)-1
                sum1 = 0
                j = idx_hopping_ion
                for k in np.arange(0,3):
                    #sum1 += np.power(v[3*j+k][i],2)
                    sum1 += np.power(np.absolute(v[3*j+k][i]),2)
                counter[int(indx)] += sum1
        # Calculate the partial PBC for the hopping O
        sum1 = 0
        sum2 = 0
        for i in np.arange(0,np.size(counter),1):
            if (w[i] > self.fmin):
                sum1 += w[i] * counter[i]
                sum2 += counter[i]
        hopping_O_PBC = sum1/sum2
        return hopping_O_PBC

    def all_atom_PBC(self, compound_name):
        N_qpoints = len(self.data['phonon'])
        N_modes_each_qpoints = len(self.data['phonon'][0]['band']) # Think of this as Nband = Nmodes
        #if (N_modes_each_qpoints != Nmodes):
        #    print('Waaaaaaarnning!!!!!!!!!!!!')
        Natoms = int(N_modes_each_qpoints/3)
        # Let's put all the frequencies in a big vector
        freq = np.zeros([N_qpoints*N_modes_each_qpoints]) # Initiating the freq array
        countfreq = 0
        for nq in np.arange(0, N_qpoints, 1):
            for nm in np.arange(0, N_modes_each_qpoints, 1):
                freq[countfreq] = self.data['phonon'][nq]['band'][nm]['frequency']
                countfreq += 1
        min_freq = min(freq)
        max_freq = max(freq)
        domega = (max_freq - min_freq)/self.Npts
        w = np.arange(1,self.Npts+1,1)*domega + min_freq + domega*.5

        # All atoms - let's go
        filename = "atom_PBC_" + compound_name + ".txt"
        f_atom_PBC = open(filename, "w")
        for atom_idx in np.arange(0, Natoms, 1):
            counter = np.zeros(self.Npts) # Histogram bins
            for nq in np.arange(0, N_qpoints, 1):
              # Let's prepare the freq vector and ev matrix for each q-point first 
                freq = np.zeros([N_modes_each_qpoints]) # Initiating the freq array
                v = np.zeros([N_modes_each_qpoints, N_modes_each_qpoints], dtype=complex) # Initiating the ev array
                for nm in np.arange(0, N_modes_each_qpoints, 1): # This loop is not needed. It will be done inside each q-point calculation.
                    freq[nm] = self.data['phonon'][nq]['band'][nm]['frequency']
                    for na in np.arange(0, Natoms, 1):
                        v.real[na*3+0, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][0][0]
                        v.real[na*3+1, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][1][0]
                        v.real[na*3+2, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][2][0]
                        v.imag[na*3+0, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][0][1]
                        v.imag[na*3+1, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][1][1]
                        v.imag[na*3+2, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][2][1]
  
            for i in np.arange(0,np.size(freq),1):
                if (np.floor((freq[i]-min_freq)/domega)<self.Npts):
                    indx = np.floor((freq[i]-min_freq)/domega)
                else:
                    indx = np.floor((freq[i]-min_freq)/domega)-1
                sum1 = 0
                j = atom_idx
                for k in np.arange(0,3):
                    #sum1 += np.power(v[3*j+k][i],2)
                    sum1 += np.power(np.absolute(v[3*j+k][i]),2)
                counter[int(indx)] += sum1
            sum1 = 0
            sum2 = 0
            for i in np.arange(0,np.size(counter),1):
                if (w[i] > self.fmin):
                    sum1 += w[i] * counter[i]
                    sum2 += counter[i]
            atom_PBC = sum1/sum2
            f_atom_PBC.write("%d %f\n" %(atom_idx+1, atom_PBC))

        f_atom_PBC.close()


    def partialDOS1(self, val1, val2): # Function to calculate partial phonon BCs
        spcs_name = val1
        AB_first = val2
        myatomNS = self.atomNs.copy()
        # Calculate the correct range for performing the summation
    
        counter = 0
        for spname in self.atomsp:
            if spname == spcs_name:
                break
            counter += 1
        begidx = 0
        endidx = 0
        for i in np.arange(0,counter+1,1):
            begidx = endidx
            endidx += myatomNS[i]
        # begidx = self.atomNs[0]-1
        # endidx = begidx + self.atomNs[1]
        #print(self.atomNs)
        print('%d %d' %(begidx, endidx))

        N_qpoints = len(self.data['phonon'])
        N_modes_each_qpoints = len(self.data['phonon'][0]['band']) # Think of this as Nband = Nmodes
        #if (N_modes_each_qpoints != Nmodes):
        #    print('Waaaaaaarnning!!!!!!!!!!!!')
        Natoms = int(N_modes_each_qpoints/3)
        #print(Natoms)
        counter = np.zeros(self.Npts) # Histogram bins
        # Let's put all the frequencies in a big vector
        freq = np.zeros([N_qpoints*N_modes_each_qpoints]) # Initiating the freq array
        countfreq = 0
        for nq in np.arange(0, N_qpoints, 1):
            for nm in np.arange(0, N_modes_each_qpoints, 1):
                freq[countfreq] = self.data['phonon'][nq]['band'][nm]['frequency']
                countfreq += 1
        min_freq = min(freq)
        max_freq = max(freq)
        domega = (max_freq - min_freq)/self.Npts
        w = np.arange(1,self.Npts+1,1)*domega + min_freq + domega*.5

        # Calculate the total DOS
        for nq in np.arange(0, N_qpoints, 1):
            # Let's prepare the freq vector and ev matrix for each q-point first 
            freq = np.zeros([N_modes_each_qpoints]) # Initiating the freq array
            v = np.zeros([N_modes_each_qpoints, N_modes_each_qpoints], dtype=complex) # Initiating the ev array
            for nm in np.arange(0, N_modes_each_qpoints, 1): # This loop is not needed. It will be done inside each q-point calculation.
                freq[nm] = self.data['phonon'][nq]['band'][nm]['frequency']
                for na in np.arange(0, Natoms, 1):
                    v.real[na*3+0, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][0][0]
                    v.real[na*3+1, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][1][0]
                    v.real[na*3+2, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][2][0]
                    v.imag[na*3+0, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][0][1]
                    v.imag[na*3+1, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][1][1]
                    v.imag[na*3+2, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][2][1]

            # Calculate the partial DOS - it is the summation of the effect coming from different q-points
            for i in np.arange(0,np.size(freq),1):
                if (np.floor((freq[i]-min_freq)/domega)<self.Npts):
                    indx = np.floor((freq[i]-min_freq)/domega)
                else:
                    indx = np.floor((freq[i]-min_freq)/domega)-1
                sum1 = 0
                for j in np.arange(begidx,endidx,1):
                    for k in np.arange(0,3):
                        #sum1 += np.power(v[3*j+k][i],2)
                        sum1 += np.power(np.absolute(v[3*j+k][i]),2)
                counter[int(indx)] += sum1
        # Calculate the partial PBC
        sum1 = 0
        sum2 = 0
        for i in np.arange(0,np.size(counter),1):
            if (w[i] > self.fmin):
                sum1 += w[i] * counter[i]
                sum2 += counter[i]
        partialPBC = sum1/sum2
        return w, counter, partialPBC
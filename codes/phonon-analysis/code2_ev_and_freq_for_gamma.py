# KG 7/5/2021

import numpy as np
import scipy
from scipy import linalg as LA

class ev_freq_info_for_gamma:
    def __init__(self,val1):
        self.mass=val1

    def check_symmetric(a, tol=1e-8):
        return np.all(np.abs(a-a.transpose()) < tol)

    def perform_LD(self,filename): # Function to calculate frequencies and eigenvectors
        # A few unit conversions
        AMU = 1.6605402e-27  # [kg]
        EV = 1.60217733e-19  # [J]
        Angstrom = 1.0e-10   # [m]
        THz = 1.0e12 # [Hz]
        factor = EV/(Angstrom**2)/np.sqrt(AMU**2) # the denominator looks like this b/c of the sqrt(m1*m2) term
        with open(f"{filename}FORCE_CONSTANTS", "r") as f:   ####
            countl = 0
            lines = f.readlines()
            templ = lines[countl].split()
            Natoms = int(templ[0])
            D = np.zeros((Natoms*3, Natoms*3))
            # print(Natoms)
            for i in range(0,Natoms,1):
                for j in range(0,Natoms,1):
                    countl = countl + 1
                    templ = lines[countl].split()
                    m = int(templ[0])
                    n = int(templ[1])
                    for k in range(0,3,1):
                        countl = countl + 1
                        templ = lines[countl].split()
                        x = float(templ[0])
                        y = float(templ[1])
                        z = float(templ[2])
                        # print(len(self.mass))
                        D[(m-1)*3+k][(n-1)*3+0] = x/np.sqrt(self.mass[m-1]*self.mass[n-1])*factor
                        D[(m-1)*3+k][(n-1)*3+1] = y/np.sqrt(self.mass[m-1]*self.mass[n-1])*factor
                        D[(m-1)*3+k][(n-1)*3+2] = z/np.sqrt(self.mass[m-1]*self.mass[n-1])*factor
        f.close()
        # Decompose the matrix to get the frequencies and eigenvectors
        #w2, self.v = np.linalg.eig(D)
        w2, self.v = LA.eigh(D) # To make sure that the symmetric nature of the D matrix is utilized
                                # https://stackoverflow.com/questions/8765310/scipy-linalg-eig-return-complex-eigenvalues-for-covariance-matrix
        self.w = w2
        countn = 0
        for x in w2:
            if x<0:
                self.w[countn] = -np.sqrt(x.astype(complex)).imag/(2*np.pi*THz)
            else:
                self.w[countn] = np.sqrt(x.astype(complex)).real/(2*np.pi*THz)
            countn += 1
        indx = np.argsort(self.w, axis=0)
        self.w = self.w[indx]
        self.v = self.v[:,indx]

        # Write the freq.txt and ev.txt files
        fw = open("freq.txt", "w")
        for i in range(0,3*Natoms,1):
            fw.write("%f " % (self.w[i]))
            fw.write("\n")

        fev = open("ev.txt", "w")
        for i in range(0,3*Natoms,1):
            for j in range(0,3*Natoms,1):
                fev.write("%f " % (self.v[i][j]))
            fev.write("\n")
        fw.close()
        fev.close()

        # Let's calculate the percentage of negative frequencies
        count_neg_freq = 0
        for i in np.arange(0, len(self.w), 1):
            if self.w[i] < 0:
                count_neg_freq += 1
        self.prcnt_neg_freq = count_neg_freq/len(self.w)*100

class phononBC:
    def __init__(self,val1, val2, val3, val4, val5):
        self.Npts = val1
        self.frequency = val2
        self.ev = val3
        self.atomNs = val4
        self.atomsp = val5

    def TotalDOS(self): # Function to calculate Total phonon BCs
        # Calculate the total DOS
        counter = np.zeros(self.Npts)
        freq = self.frequency.copy()
        min_freq = min(freq)
        max_freq = max(freq)
        domega = (max_freq - min_freq)/self.Npts
        w = np.arange(1,self.Npts+1,1)*domega + min_freq + domega*.5
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
        for i in np.arange(0,counter+1,1):
            begidx = endidx
            endidx += myatomNS[i]
        #print('%d %d' %(begidx, endidx))

        # Calculate the partial DOS
        counter = np.zeros(self.Npts)
        freq = self.frequency.copy()
        v = self.ev.copy()
        min_freq = min(freq)
        max_freq = max(freq)
        domega = (max_freq - min_freq)/self.Npts
        w = np.arange(1,self.Npts+1,1)*domega + min_freq + domega*.5
        for i in np.arange(0,np.size(freq),1):
            if (np.floor((freq[i]-min_freq)/domega)<self.Npts):
                indx = np.floor((freq[i]-min_freq)/domega)
            else:
                indx = np.floor((freq[i]-min_freq)/domega)-1
            sum1 = 0
            for j in np.arange(begidx,endidx,1):
                for k in np.arange(0,3):
                    sum1 += np.power(v[3*j+k][i],2)
            counter[int(indx)] += sum1
        # Calculate the partial PBC
        sum1 = 0
        sum2 = 0
        for i in np.arange(0,np.size(counter),1):
            sum1 += w[i] * counter[i]
            sum2 += counter[i]
        partialPBC = sum1/sum2
        return w, counter, partialPBC

    def ave_freq_hopping_O(self, val1): # Function to calculate the average frequency of the hopping ion

        idx_hopping_ion = val1
        counter = np.zeros(self.Npts)
        freq = self.frequency.copy()
        v = self.ev.copy()
        min_freq = min(freq)
        max_freq = max(freq)
        domega = (max_freq - min_freq)/self.Npts
        w = np.arange(1,self.Npts+1,1)*domega + min_freq + domega*.5
        for i in np.arange(0,np.size(freq),1):
            if (np.floor((freq[i]-min_freq)/domega)<self.Npts):
                indx = np.floor((freq[i]-min_freq)/domega)
            else:
                indx = np.floor((freq[i]-min_freq)/domega)-1
            sum1 = 0
            j = idx_hopping_ion
            for k in np.arange(0,3):
                sum1 += np.power(v[3*j+k][i],2)
            counter[int(indx)] += sum1
        # Calculate the partial PBC for the hopping O
        sum1 = 0
        sum2 = 0
        for i in np.arange(0,np.size(counter),1):
            sum1 += w[i] * counter[i]
            sum2 += counter[i]
        hopping_O_PBC = sum1/sum2
        return hopping_O_PBC


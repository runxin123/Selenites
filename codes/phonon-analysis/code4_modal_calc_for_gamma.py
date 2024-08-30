# KG 7/6/2021

import numpy as np

class modal_calc_for_gamma:
    def __init__(self,val1, val2, val3, val4, val5):
        self.dx_cartesian = val1 # Displacement field in Cartesian coordinates
        self.w = val2 # Phonon frequencies
        self.v = val3 # Phonon eigenvectors
        self.mass = val4 # Atomic masses
        self.Natoms = val5 # Number of atoms

    def perform_modal(self):
        # A few unit conversions
        AMU = 1.6605402e-27 # [kg]
        EV = 1.60217733e-19 # [J]
        Angstrom = 1.0e-10  # [m]
        THz = 1.0e12 # [Hz]

        Nmodes = 3 * self.Natoms # Number os modes
        Qn = np.zeros(Nmodes) # Qn stores modal displacements
        self.En = np.zeros(Nmodes) # En stores modal excitations (energies)
        for n in np.arange(0,Nmodes,1): # Looping over number of modes
            sum1 = 0
            for i in np.arange(0,self.Natoms,1): # Looping over number of atoms
                sum1 += np.sqrt(self.mass[i]*AMU)*np.dot(self.v[i*3+0:i*3+3,n],(self.dx_cartesian[i,:]*Angstrom))
            Qn[n] = sum1
            self.En[n] = 0.5*np.power(self.w[n]*THz*2*np.pi,2)*np.power(Qn[n],2)/EV # in eV
        '''fEn = open("En.txt", "w")
        for n in range(0,Nmodes,1):
            fEn.write("%f " % (self.En[n]))
            fEn.write("\n")
        fEn.close()'''

class En_hist_n_center_for_gamma:
    def __init__(self,val1, val2, val3):
        self.freq = val1 # Phonon frequencies
        self.En = val2 # Number of bins
        self.Npts = val3 # Number of bins
    def hist_for_En_n_center(self):
        # Calculate En histogram
        counter = np.zeros(self.Npts)
        freq = self.freq.copy()
        En = self.En.copy()
        min_freq = min(freq)
        max_freq = max(freq)
        domega = (max_freq - min_freq)/self.Npts
        w = np.arange(1,self.Npts+1,1)*domega + min_freq + domega*.5
        for i in np.arange(1,np.size(freq),1):
            if (np.floor((freq[i]-min_freq)/domega)<self.Npts):
                indx = np.floor((freq[i]-min_freq)/domega)
            else:
                indx = np.floor((freq[i]-min_freq)/domega)-1
            counter[int(indx)] += En[i];
        # Calculate En center
        sum1 = 0
        sum2 = 0
        for i in np.arange(0,np.size(counter),1):
            sum1 += w[i] * counter[i]
            sum2 += counter[i]
        w_En_center = sum1/sum2
        return w, counter, w_En_center

    def freq_max_En(self):
        # the frequency associated with maximum En
        return(self.freq[np.argmax(self.En)])


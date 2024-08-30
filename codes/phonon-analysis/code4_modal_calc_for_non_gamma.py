# KG 7/6/2021

import numpy as np

class modal_calc:
    def __init__(self,val1, val2, val3, val4):
        self.dx_cartesian = val1 # Displacement field in Cartesian coordinates
        self.data = val2 # Phonon information
        self.mass = val3 # Atomic masses
        self.Natoms = val4 # Number of atoms

    def perform_modal(self):
        # A few unit conversions
        AMU = 1.6605402e-27 # [kg]
        EV = 1.60217733e-19 # [J]
        Angstrom = 1.0e-10  # [m]
        THz = 1.0e12 # [Hz]

        rlattice_vectors = np.asarray(self.data['reciprocal_lattice']) # Reciprocqal lattice vectors
        lattice_vectors = np.asarray(self.data['lattice']) # lattice vectors
        #pos = np.asarray(self.data['points'][0]['coordinates'][0]))

        Nmodes = 3 * self.Natoms # Number os modes
        N_qpoints = len(self.data['phonon'])
        N_modes_each_qpoints = len(self.data['phonon'][0]['band']) # Think of this as Nband = Nmodes

        Qn = np.zeros(Nmodes, dtype=complex) # Qn stores modal displacements
        self.En = np.zeros([Nmodes, N_qpoints]) # En stores modal excitations (energies)

        if (N_modes_each_qpoints != Nmodes):
            print('Waaaaaaarnning1!!!!!!!!!!!!')
        for nq in np.arange(0, N_qpoints, 1):
            # Let's prepare the freq vector and ev matrix for each q-point first 
            w = np.zeros([N_modes_each_qpoints]) # Initiating the freq array
            v = np.zeros([N_modes_each_qpoints, N_modes_each_qpoints], dtype=complex) # Initiating the ev array
            for nm in np.arange(0, N_modes_each_qpoints, 1): # This loop is not needed. It will be done inside each q-point calculation.
                w[nm] = self.data['phonon'][nq]['band'][nm]['frequency']
                for na in np.arange(0, self.Natoms, 1):
                    v.real[na*3+0, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][0][0]
                    v.real[na*3+1, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][1][0]
                    v.real[na*3+2, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][2][0]
                    v.imag[na*3+0, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][0][1]
                    v.imag[na*3+1, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][1][1]
                    v.imag[na*3+2, nm] = self.data['phonon'][nq]['band'][nm]['eigenvector'][na][2][1]

            q_position = np.asarray(self.data['phonon'][nq]['q-position'])
            wave_vector = np.matmul(q_position, rlattice_vectors)
            for n in np.arange(0,Nmodes,1): # Looping over number of modes
                sum1 = 0
                for i in np.arange(0,self.Natoms,1): # Looping over number of atoms
                    pos_direct = np.asarray(self.data['points'][0]['coordinates'])
                    pos_cartesian = np.matmul(pos_direct, lattice_vectors)
                    kdotr = np.zeros((1), dtype=complex)
                    kdotr.imag = -np.dot(wave_vector, pos_cartesian)
                    expkr = np.exp(kdotr)
                    sum1 += np.sqrt(self.mass[i]*AMU)*expkr*np.dot(np.conjugate(v[i*3+0:i*3+3,n]),self.dx_cartesian[i,:]*Angstrom)
                Qn[n] = sum1
                self.En[n][nq] = np.real(0.5*np.power(w[n]*THz*2*np.pi,2)*np.conjugate(Qn[n])*Qn[n]/EV) # in eV
        '''fEn = open("En.txt", "w")
        for n in range(0,Nmodes,1):
            fEn.write("%f " % (self.En[n]))
            fEn.write("\n")
        fEn.close()'''

class En_hist_n_center:
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


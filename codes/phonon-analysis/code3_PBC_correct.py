# KG 7/6/2021

import numpy as np
import copy

class PBC_correct:
    def __init__(self,val1, val2, val3, val4):
        self.dx_direct = val1
        self.L = val2
        self.coeff = val3
        self.Natoms = val4

    def perform_PBC_correct(self):
        L = self.L.copy()
        L = L * self.coeff
        # Calculate the unit vectors along the 1st, 2nd, and 3rd main axes.
        L0 = np.sqrt(np.sum(np.power(L[0,:], 2)))
        L1 = np.sqrt(np.sum(np.power(L[1,:], 2)))
        L2 = np.sqrt(np.sum(np.power(L[2,:], 2)))
        L0hat = L[0,:]/L0
        L1hat = L[1,:]/L1
        L2hat = L[2,:]/L2
        # Convert dx from direct to Cartesian
        self.dx_cartesian  = self.dx_direct @ L
        # Perform PBC
        for i in range(self.Natoms):
            # correct along the 1st vector
            d0 = np.dot(self.dx_cartesian[i,:],L0hat)
            if (d0 >= L0/2):
                self.dx_cartesian[i,:] -= L[0,:]
            elif (d0 < -L0/2):
                self.dx_cartesian[i,:] += L[0,:]
            # correct along the 2nd vector
            d1 = np.dot(self.dx_cartesian[i,:],L1hat)
            if (d1 >= L1/2):
                self.dx_cartesian[i,:] -= L[1,:]
            elif (d1 < -L1/2):
                self.dx_cartesian[i,:] += L[1,:]
            # correct along the 3rd vector
            d2 = np.dot(self.dx_cartesian[i,:],L2hat)
            if (d2 >= L2/2):
                self.dx_cartesian[i,:] -= L[2,:]
            elif (d2 < -L2/2):
                self.dx_cartesian[i,:] += L[2,:]


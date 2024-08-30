# KG 7/5/2021

import os
import numpy as np
import pymatgen.core as pmg

class POSCARs_info_NEB:
    def __init__(self, val1, val2, val3, val4, val5):
        self.idx1 = val1
        self.idx2 = val2
        self.idx3 = val3
        self.idx4 = val4
        self.path_NEB = val5

    def read_masses_and_xyzs(self):
        filename = os.path.join(self.path_NEB, "0" + str(self.idx1),"POSCAR")
        if os.path.exists(filename) == False:
            filename = os.path.join(self.path_NEB, "0" + str(self.idx1),"CONTCAR")
        self.L = np.array([[0, 0, 0]])
        self.x1 = np.array([[0, 0, 0]])
        # Open POSCAR 1
        with open(filename, "r") as f:
            lines = f.readlines()
            # Read the coefficient
            templ = lines[1].split() 
            self.coeff = float(templ[0])
            # Read the cell vectors
            for i in range(2,5,1):
                templ = lines[i].split() 
                self.L = np.append(self.L, [[float(templ[0]),float(templ[1]),float(templ[2])]],axis = 0)
            self.L = np.delete(self.L, 0, 0)
            # Read the atom species
            templ = lines[5].split() 
            atomsp = []
            Nspecies = len(templ)
            mass_sp = []
            for i in range(0,Nspecies,1):
                atomsp.append(templ[i])
                mass_sp.append(float(pmg.Element(atomsp[i]).atomic_mass))
            self.atomsp_passed = atomsp
            # Read the atom numbers for each species
            templ = lines[6].split()
            atomNs = np.empty(0, dtype='int')
            for i in range(0,Nspecies,1):
                atomNs = np.append(atomNs, int(templ[i]))
            self.atomNs_passed = atomNs
            self.Natoms = np.sum(atomNs)
            # create the mass array
            self.mass = np.array([])
            for i in range(0,Nspecies,1):
                self.mass = np.append(self.mass, np.ones(atomNs[i])*mass_sp[i])
            # Read direct or cartesian
            #templ = lines[7].split()
            #direct_or_cartesian = 
            # Read image 1 coordinates
            for i in range(8,8+self.Natoms,1):
                templ = lines[i].split() 
                self.x1 = np.append(self.x1, [[float(templ[0]),float(templ[1]),float(templ[2])]],axis = 0)
            self.x1 = np.delete(self.x1, 0, 0)
        f.close()

        # Open POSCAR 2
        # Just read the coordinates
        filename = os.path.join(self.path_NEB, "0" + str(self.idx2), "POSCAR")
        if os.path.exists(filename) == False:
            filename = os.path.join(self.path_NEB, "0" + str(self.idx2), "CONTCAR")
        self.x2 = np.array([[0, 0, 0]])
        # Open POSCAR 1
        with open(filename, "r") as f:
            lines = f.readlines()
            for i in range(8,8+self.Natoms,1):
                templ = lines[i].split() 
                self.x2 = np.append(self.x2, [[float(templ[0]),float(templ[1]),float(templ[2])]],axis = 0)
            self.x2 = np.delete(self.x2, 0, 0)
        f.close()

        # Open POSCAR 3 (corresponding to TS image)
        # Just read the coordinates
        filename = os.path.join(self.path_NEB, "0" + str(self.idx3), "POSCAR")
        if os.path.exists(filename) == False:
            filename = os.path.join(self.path_NEB, "0" + str(self.idx3), "CONTCAR")
        self.x3 = np.array([[0, 0, 0]])
        # Open POSCAR 1
        with open(filename, "r") as f:
            lines = f.readlines()
            for i in range(8,8+self.Natoms,1):
                templ = lines[i].split() 
                self.x3 = np.append(self.x3, [[float(templ[0]),float(templ[1]),float(templ[2])]],axis = 0)
            self.x3 = np.delete(self.x3, 0, 0)
        f.close()

        # Open POSCAR final (image 4 usually) (corresponding to final image)
        # Just read the coordinates
        filename = os.path.join(self.path_NEB, "0" + str(self.idx4), "POSCAR")
        if os.path.exists(filename) == False:
            filename = os.path.join(self.path_NEB, "0" + str(self.idx4), "CONTCAR")
        self.x4 = np.array([[0, 0, 0]])
        # Open POSCAR 1
        with open(filename, "r") as f:
            lines = f.readlines()
            for i in range(8,8+self.Natoms,1):
                templ = lines[i].split() 
                self.x4 = np.append(self.x4, [[float(templ[0]),float(templ[1]),float(templ[2])]],axis = 0)
            self.x4 = np.delete(self.x4, 0, 0)
        f.close()

class POSCARs_info_non_NEB:
    def __init__(self, val1):
        self.POSCAR_path = val1

    def read_masses_and_xyzs(self):
        filename = self.POSCAR_path + "/POSCAR"
        self.L = np.array([[0, 0, 0]])
        self.x1 = np.array([[0, 0, 0]])
        # Open POSCAR
        with open(filename, "r") as f:
            lines = f.readlines()
            # Read the coefficient
            templ = lines[1].split() 
            self.coeff = float(templ[0])
            # Read the cell vectors
            for i in range(2,5,1):
                templ = lines[i].split() 
                self.L = np.append(self.L, [[float(templ[0]),float(templ[1]),float(templ[2])]],axis = 0)
            self.L = np.delete(self.L, 0, 0)
            # Read the atom species
            templ = lines[5].split() 
            atomsp = []
            Nspecies = len(templ)
            mass_sp = []
            for i in range(0,Nspecies,1):
                atomsp.append(templ[i])
                mass_sp.append(float(pmg.Element(atomsp[i]).atomic_mass))
            self.atomsp_passed = atomsp
            # Read the atom numbers for each species
            templ = lines[6].split()
            atomNs = np.empty(0, dtype='int')
            for i in range(0,Nspecies,1):
                atomNs = np.append(atomNs, int(templ[i]))
            self.atomNs_passed = atomNs
            self.Natoms = np.sum(atomNs)
            # create the mass array
            self.mass = np.array([])
            for i in range(0,Nspecies,1):
                self.mass = np.append(self.mass, np.ones(atomNs[i])*mass_sp[i])
            # Read direct or cartesian
            #templ = lines[7].split()
            #direct_or_cartesian = 
            # Read image 1 coordinates
            for i in range(8,8+self.Natoms,1):
                templ = lines[i].split() 
                self.x1 = np.append(self.x1, [[float(templ[0]),float(templ[1]),float(templ[2])]],axis = 0)
            self.x1 = np.delete(self.x1, 0, 0)
        f.close()


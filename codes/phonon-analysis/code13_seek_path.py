# KG 11/1/2021

import numpy as np
import seekpath

def BZ_vectors_from_package(structure, with_time_reversal, recipe, threshold, symprec, angle_tolerance):
    skpth = seekpath.get_path(structure, with_time_reversal, recipe, threshold, symprec, angle_tolerance)
    #print(skpth)
    #L = skpth['primitive_lattice']
    #x = skpth['primitive_positions']
    L = skpth['conv_lattice']
    x = skpth['conv_positions']
    rl_vectors = skpth['reciprocal_primitive_lattice']
    return L, x, rl_vectors

def BZ_vectors_traditional_calc(L): # https://en.wikipedia.org/wiki/Reciprocal_lattice#Three_dimensions
    a1 = L[0,:]
    a2 = L[1,:]
    a3 = L[2,:]
    V = np.dot(a1,np.cross(a2,a3))
    rl_vectors = np.zeros([3, 3])
    #rl_vectors[0,:] = 2*np.pi/V * np.cross(a2,a3)
    #rl_vectors[1,:] = 2*np.pi/V * np.cross(a3,a1)
    #rl_vectors[2,:] = 2*np.pi/V * np.cross(a1,a2)
    #Let's multiply 2*pi later in the dynamical matrix formation
    rl_vectors[0,:] = 1./V * np.cross(a2,a3)
    rl_vectors[1,:] = 1./V * np.cross(a3,a1)
    rl_vectors[2,:] = 1./V * np.cross(a1,a2)

    return rl_vectors




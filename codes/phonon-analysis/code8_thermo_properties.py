# KG 7/19/2021

import numpy as np
from math import exp

# A few unit conversions
AMU = 1.6605402e-27  # [kg]
Angstrom = 1.0e-10   # [m]
THz = 1.0e12 # [Hz]
hbar = 6.63e-34/(2*np.pi) # J.s
kB = 1.38064852e-23 # J/K
EV = 1.60217733e-19 # [J]
NAvo = 6.0221409e+23 # Avogadro's number

def vibrational_entropy_diff(w1, w2, Natoms, T, Natoms_in_compound_formula):

    # vibrational entropy formula will break with negative frequencies
    # vibrational entropy formula here: https://phonopy.github.io/phonopy/formulation.html#entropy
    w1 = np.abs(w1)
    w2 = np.abs(w2)

    Nmodes = 3*Natoms
    S1 = 0
    S2 = 0
    for n in np.arange(0,Nmodes,1):
        w_rad1 = 2*np.pi*THz*w1[n]
        x1 = hbar*w_rad1 / (2*kB*T)

        w_rad2 = 2*np.pi*THz*w2[n]
        x2 = hbar*w_rad2 / (2*kB*T)

        S1 += 1./(2*T) * hbar*w_rad1 * 1./np.tanh(x1) - kB * np.log(2*np.sinh(x1))
        S2 += 1./(2*T) * hbar*w_rad2 * 1./np.tanh(x2) - kB * np.log(2*np.sinh(x2))

    diff_vib_S = S2 - S1 # J/K
    diff_vib_S_percent_change = diff_vib_S/S1*100. # percent
    #diff_vib_S /= EV # eV/K
    diff_vib_S *= 1./(Natoms/(Natoms_in_compound_formula * NAvo)) # J/mol-K # denominator is the number of atoms that are present in one mole of compound
    return diff_vib_S, diff_vib_S_percent_change

# vibrational energy
# add another column based on zero frequency

def vibrational_energy_diff(w1, w2, Natoms, T):

    # vibrational energy formula will be weird with negative frequencies
    # vibrational energy formula here: https://phonopy.github.io/phonopy/formulation.html#harmonic-phonon-energy
    w1 = np.abs(w1)
    w2 = np.abs(w2)

    Nmodes = 3*Natoms
    E1 = 0
    E2 = 0
    for n in np.arange(0,Nmodes,1):
        w_rad1 = 2*np.pi*THz*w1[n]
        x1 = hbar*w_rad1 / (kB*T)

        w_rad2 = 2*np.pi*THz*w2[n]
        x2 = hbar*w_rad2 / (kB*T)

        E1 += hbar*w_rad1 * (1./2 + 1./(exp(x1)-1))
        E2 += hbar*w_rad2 * (1./2 + 1./(exp(x2)-1))

    E1 /= EV # eV
    E2 /= EV # eV

    diff_vib_E = E2 - E1 # eV
    diff_vib_E_percent_change = diff_vib_E/E1*100. # percent

    return diff_vib_E, diff_vib_E_percent_change


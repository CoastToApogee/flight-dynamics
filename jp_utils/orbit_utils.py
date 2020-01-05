"""
orbit_utils
Utilties related to orbits. Calculating various parameters, converting between
different coordinates, etc.
"""

import numpy as np
from jp_utils.general import *
from jp_utils.constants.earth import *

def calc_energy(state_vector):
    """ Calculate energy from a state vector.
    """

    (pos, vel) = np.split(state_vector, 2)
    return np.sum(vel**2) / 2.0 - MU/np.sum(pos**2)**0.5

def calc_h(state_vector):
    """ Calculate angular momentum vector from a state vector.
    """

    (pos, vel) = np.split(state_vector, 2)
    return np.cross(pos, vel)

def calc_ecc(state_vector):
    """ Calculate eccentricity vector from a state vector.
    """

    (pos, vel) = np.split(state_vector, 2)

    h = calc_h(state_vector)

    return np.cross(vel, h)/MU - pos / np.sum(pos**2)**0.5

def cartesian_to_keplerian(state):
    """ Convert a cartesian state vector into a keplerian state vector.

    Note that the keplerian vector is a dictionary keyed on Keplerian
    parameter names.
    """

    (pos, vel) = np.split(state, 2)

    h = np.cross(pos, vel)

    ecc_vec = np.cross(vel, h) / MU - (pos / vec_mag(pos))

    n = np.cross([0,0,1], h)

    ecc = vec_mag(ecc_vec)

    inc = np.arccos(h[2]/vec_mag(h))

    ran = np.arccos(n[0]/vec_mag(n))

    aop = np.arccos(np.dot(n,ecc_vec) / (vec_mag(n) * ecc))

    sma = 1 / (2/vec_mag(pos) - vec_mag(vel)**2/MU)

    tam = np.arccos(np.dot(ecc_vec, pos) / (ecc * vec_mag(pos)))

    return {
        'inc' : inc,
        'ecc' : ecc,
        'ran' : ran,
        'aop' : aop,
        'sma' : sma,
        'tam' : tam,
    }

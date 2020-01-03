"""
orbit_utils
Utilties related to orbits. Calculating various parameters, converting between
different coordinates, etc.
"""

import numpy as np
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

    ecc_vec = np.cross(vel, h) / MU - (pos / norm(pos))
    ecc = norm(ecc_vec)

    n = np.cross([0,0,1], h)

    inc = np.arccos(h[2]/norm(h))

    ran = np.arccos(n[0]/norm(n))

    aop = np.arccos(np.dot(n,ecc_vec) / (norm(n) * ecc))

    sma = 1 / (2/norm(pos) - norm(vel)**2/MU)

    return {
        'inc' : inc,
        'ecc' : ecc,
        'ran' : ran,
        'aop' : aop,
        'sma' : sma
    }

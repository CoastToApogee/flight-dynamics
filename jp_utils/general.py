"""
general
A place to park stuff until I think of a better place to put it.
"""

import numpy as np
from jp_utils.constants.earth import *

def accel_point_mass(t, state):
    (pos, vel) = np.split(state, 2)
    return np.array([ vel, -(MU / np.sum(pos**2)**1.5) * pos ]).flatten()

def vec_mag(vec):
    return np.sqrt( (vec**2).sum() )

def accel_j2(t, state):
    # from curtis 12.30
    (pos, vel) = np.split(state, 2)
    r = vec_mag(pos)
    (x,y,z) = pos

    scalar = 1.5 * (J2 * MU * RADIUS**2)/r**4
    z2_over_r2 = 5*z**2/r**2
    vector = np.array([
        x/r * (z2_over_r2 - 1),
        y/r * (z2_over_r2 - 1),
        z/r * (z2_over_r2 - 3),
    ])

    acc = scalar * vector
    return np.array([0,0,0, acc[0], acc[1], acc[2]])

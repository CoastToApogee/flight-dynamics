import numpy as np

def accel_point_mass(t, state):
    (pos, vel) = np.split(state, 2)
    return np.array([ vel, -(MU / np.sum(pos**2)**1.5) * pos ]).flatten()

def accel_j2(t, state):
    # from curtis 12.30
    (pos, vel) = np.split(state, 2)
    r = np.linalg.norm(pos)
    (x,y,z) = pos

    scalar = 1.5 * (J2 * MU * EARTH_RADIUS**2)/r**4
    z2_over_r2 = 5*z**2/r**2
    vector = np.array([
        x/r * (z2_over_r2 - 1),
        y/r * (z2_over_r2 - 1),
        z/r * (z2_over_r2 - 3),
    ])
    return np.array([0,0,0], scalar * vector)

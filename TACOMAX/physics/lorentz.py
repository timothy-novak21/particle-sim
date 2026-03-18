import numpy as np


def accel(t,state,field,CONST):
    """
    Compute the acceleration on a charged particle as a result of Lorentz force

    Parameters
    ----------
    t : scalar
        Current time in seconds
        Note: t is required by solve_ivp to integrate but is unused because acceleration is independent of the time

    state : array, shape (6,)
        Particle state vector, first 3 entries are cartesian position in meters, next 3 entires are cartesian velocity components in meters per second

    field : object
        Magnetic field object with a vector(pos) method
        Compatible with Mirror and Tokamak classes

    CONST : named tuple, shape (2,)
        Named tuple containing the physical constants q (proton charge) in Coulombs and m (proton mass) in kilograms

    Returns
    -------
    delta_state : array, shape (6,)
        Delta state function of the particle subject to Lorentz forces
    """

    pos = state[0:3] # position component of state vector [m]
    vel = state[3:6] # velocity component of state vector [m/s]

    B_pos = field.vector(pos) # B field at particle position [T]
    a = (CONST.q / CONST.m) * np.cross(vel, B_pos) # a = F/m = (q/m)*(E + v X B) [m/s^2]

    delta_state = np.concatenate([vel, # change in position is the velocity of the particle [m/s]
                                  a])  # change in velocity is the acceleration of the particle [m/s]

    return delta_state
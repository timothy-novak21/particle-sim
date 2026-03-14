import numpy as np

def check_mu_conservation(state,mirror,CONST):
    """
    Check percent change in magnetic moment to determine accuracy of simulation

    Parameters
    ----------
    state : array, shape (6,)
        Particle state vector, first 3 entries are cartesian position in meters, next 3 entires are cartesian velocity components in meters per second

    mirror : object
        Magnetic mirror object belonging to the Mirror class stored in B_field.py

    CONST : named tuple, shape (2,)
        Named tuple containing the physical constants q (proton charge) in Coulombs and m (proton mass) in kilograms

    Returns
    -------
    drift : scalar
        Percent change in magnetic moment to validate fidelity of simulation, unitless
    """

    mus = []
    for i in range(state.y.shape[1]):
        pos = state.y[0:3, i]
        vel = state.y[3:6, i]
        B_vec = mirror.vector(pos)
        B_mag = np.linalg.norm(B_vec)
        B_hat = B_vec / B_mag

        v_par_vec = np.dot(vel,B_hat) * B_hat
        v_perp_vec = vel - v_par_vec
        v_perp_sq = np.dot(v_perp_vec, v_perp_vec)

        mu = 0.5 * CONST.m * v_perp_sq / B_mag
        mus.append(mu)

    drift = (np.max(mus) - np.min(mus)) / mus[0]

    return drift
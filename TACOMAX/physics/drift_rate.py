import numpy as np


def check_drift_rate(state,toroid,v0,theta,CONST):
    """
    Measure simulated vertical drift rate and compare to analytic drift rate prediction

    Parameters
    ----------
    state : OdeResult
        Solution object from solve_ivp
    
    toroid : Tokamak object
        Magnetic field object belonging to the Tokamak class stored in tokamak.py

    v0 : scalar
        Initial particle speed in meters per second

    theta : scalar
        Particle pitch angle in radians

    CONST : named tuple, shape (2,)
        Named tuple containing the physical constants q (proton charge) in Coulombs and m (proton mass) in kilograms
        
    Returns
    -------
    drift_sim : scalar
        Simulated vertical drift rate in meters per second

    drift_analytic : scalar
        Analytic vertical drift rate prediction in meters per second

    err : scalar
        Percent error between drift rate simulated and analytic prediction, unitless
    """

    # Calculate analytic grad-B drift
    v_perp = v0 * np.sin(theta) # perpendicular velocity [m/s]
    v_par = v0 * np.cos(theta) # parallel velocity [m/s]
    omega_c = CONST.q * toroid.B0 / CONST.m # cyclotron frequency [1/s]
    r_L = v_perp / omega_c # Larmor radius [m]

    v_gradB = CONST.m * v_perp**2 / (2 * CONST.q * toroid.B0 * toroid.R0) * (1 + r_L**2 / toroid.R0**2) # grad B drift term [m/s]
    v_curvature = CONST.m * v_par**2 / (CONST.q * toroid.B0 * toroid.R0) # curvature drift term [m/s]
    drift_analytic = v_gradB + v_curvature # total vertical drift [m/s]

    # Find simulated grad-B drift
    t_sim = state.t
    z_sim = state.y[2]

    # Drift rate is the slope of a linear fit line through z position over time
    fit = np.polyfit(t_sim,z_sim,deg=1)
    drift_sim = fit[0]

    # Percent error in sim drift
    err = abs(drift_analytic - abs(drift_sim)) / drift_analytic * 100

    return drift_sim, drift_analytic, err
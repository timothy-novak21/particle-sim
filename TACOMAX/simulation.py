import numpy as np
import scipy.integrate as integrate

from fields.magnetic_mirror import Mirror
from fields.tokamak import Tokamak
from physics.lorentz import accel
from analysis.plotting import plot_pitch_sweep, plot_crit_curve, plot_drift, plot_toroid_poloid_comp, plot_q_sweep
from physics.adiabaticity import invariance_check
from physics.mu_conservation import check_mu_conservation
from physics.drift_rate import check_drift_rate


def sim_single_mirror(mirror,theta,v0,mode,CONST):
    """
    Run a single particle magnetic mirror simulation

    Parameters
    ----------
    mirror : Mirror object
        Magnetic mirror object belonging to the Mirror class stored in B_field.py

    theta : scalar
        Particle pitch angle in radians

    v0 : scalar
        Initial particle speed in meters per second

    mode : string
        "sweep" | "mu" | "single"
        String containing which event set to run ODE solver with

    CONST : named tuple, shape (2,)
        Named tuple containing the physical constants q (proton charge) in Coulombs and m (proton mass) in kilograms
        
    Returns
    -------
    state : OdeResult
        Solution object from solve_ivp. state.t_events[0] == 1 indicates particle escape
        state.t_events[1] == 1 indicates particle was reflected and remained in the mirror
    """

    # Initial conditions
    v0_perp = v0 * np.sin(theta) # initial particle velocity perpendicular to B field [m/s]
    v0_par = v0 * np.cos(theta) # initial particle velocity parallel to B field [m/s]

    # Package initial state vector
    state0 = [0.0, 0.0, 0.0, # initial position [m]
              v0_perp, 0.0, v0_par] # initial velocity [m/s]
    
    # Integration time conditions computed based on particle velocity and mirror length
    T_ref_est = 4 * mirror.length / v0 # estimate reflection period [s] 
    t_span = [0,4 * T_ref_est] # integration time bounds [s]

    omega_c = CONST.q * mirror.Bmax / CONST.m # cyclotron angular frequency [rad/s]
    T_cyclotron = 2 * np.pi / omega_c # cyclotron period [s]
    max_step = T_cyclotron / 200 # 200 steps per gyration for max timestep [s]

    # Define event to terminate solver if particle escapes mirror
    def escaped(t,state,mirror,CONST):
        return abs(state[2]) - mirror.length # particle crosses when z is at either mirror boundary

    escaped.terminal = True # stop integration upon escape
    escaped.direction = 1 # stop integration when moving outward

    # Define event to terminate solver if particle is reflected (v_par = 0)
    def reflected(t,state,mirror,CONST):
        return state[5] # particles reflects when v_par changes sign
    
    reflected.terminal = True # stop integration if particle v_par changes sign
    reflected.direction = -1 # stop integration when v_par crosses zero going negative
    
    # Solve ODE with events set based on mode
    if mode == "sweep":
        # sweep mode runs with both escaped and reflected termination events
        state = integrate.solve_ivp(accel, t_span, state0,
                                    method="RK45", max_step=max_step,
                                    rtol=1e-10, atol=1e-12,
                                    args=(mirror,CONST),
                                    events=[escaped,reflected])
    elif mode == "mu":
        # mu mode runs with only the escaped termination event
        T_bounce_mu = 4 * mirror.length / v0_par
        t_span_mu = [0, 5 * T_bounce_mu] # 5 period timespan for mu conservation check
        
        state = integrate.solve_ivp(accel, t_span_mu, state0,
                                    method="RK45", max_step=max_step,
                                    rtol=1e-10, atol=1e-12,
                                    args=(mirror,CONST),
                                    events=[escaped])
    elif mode == "single":
        # single mode runs with no termination events
        state = integrate.solve_ivp(accel, t_span, state0,
                                    method="RK45", max_step=max_step,
                                    rtol=1e-10, atol=1e-12,
                                    args=(mirror,CONST))
    else:
        raise ValueError(f"Invalid mode '{mode}'. Expected 'sweep', 'mu', or 'single'.")

    return state


def sim_mirror_sweep(mirror,v0,N,CONST):
    """
    Run a series of simulations sweeping through a collection of pitch angles between 1 and 89 degrees
    Calculate and return a simulated critical pitch angle
    
    Parameters
    ----------
    mirror : Mirror object
        Magnetic mirror object belonging to the Mirror class stored in B_field.py

    v0 : scalar
        Initial particle speed in meters per second

    N : scalar
        Number of discrete pitch angles to sweep through between 1 and 89 degrees, unitless

    CONST : named tuple, shape (2,)
        Named tuple containing the physical constants q (proton charge) in Coulombs and m (proton mass) in kilograms
        
    Returns
    -------
    theta_c_sim : scalar
        Simulated critical pitch angle for a given mirror ratio in degrees
        Particles with pitch angles beneath this value will escape the mirror
    
    err_theta_c : scalar
        Percent error in simulated critical pitch angle vs analytical critical pitch angle, unitless
    """

    theta = np.linspace(np.radians(1), np.radians(89), N) # particle pitch angle [rad]
    esc = np.zeros(len(theta)) # initialize an array to track which pitch angles escape

    # Run simulations
    mode = "sweep"
    for i in range(0,len(theta)):
        print(f"Running simulation {i+1}/{N}: theta = {np.degrees(theta[i]):.1f} deg")
        state = sim_single_mirror(mirror,theta[i],v0,mode,CONST)
        esc[i] = len(state.t_events[0]) > 0 # t_events[0] == 1 when particle escapes mirror
        
    # Calc simulated critical pitch angle to compare to analytical
    if np.any(esc == 0) and np.any(esc == 1):
        high_esc = theta[esc == 1][-1] # find highest angle where particle escapes [rad]
        low_ref = theta[esc == 0][0] # find lowest angle where particle is reflected [rad]
        theta_c_sim = np.degrees((high_esc + low_ref) / 2) # [deg]
    else:
        print("Warning: no clean transition found in sweep")
        theta_c_sim = None

    # Calc analytical critical pitch angle
    theta_c_analytical = mirror.theta_c # [deg]
    err_theta_c = np.abs(theta_c_analytical - theta_c_sim) / theta_c_analytical * 100 # percent error in sim critical pitch angle [-]

    [epsilon,adiabatic] = invariance_check(v0,mirror,CONST) # adiabaticity condition for given mirror and v0

    # Plot particle end state (escaped or reflected) as a function of pitch angle
    plot_pitch_sweep(theta,esc,err_theta_c,epsilon,adiabatic,mirror)
    
    return [theta_c_sim,err_theta_c]


def sim_rm_sweep(Rm,Bmax,v0,N,CONST):
    """
    Run through pitch sweeps for each mirror ratio and plot critical pitch angle curve
    
    Parameters
    ----------
    Rm : array-like, shape (L,)
        List of mirror ratios, L terms long, unitless

    Bmax : scalar
        Maximum magnetic field strength in Tesla

    v0 : scalar
        Initial particle speed in meters per second

    N : scalar
        Number of discrete pitch angles to sweep through between 1 and 89 degrees, unitless

    CONST : named tuple, shape (2,)
        Named tuple containing the physical constants q (proton charge) in Coulombs and m (proton mass) in kilograms
        
    Returns
    -------
    none
    """

    # Set up list of Mirror objects
    mirror_list = []
    for i in range(0,len(Rm)):
        length_i = Rm[i] * 0.5
        temp_mirror = Mirror(length_i,Rm[i],Bmax) # create mirror object for current Rm
        mirror_list.append(temp_mirror) # add to mirror list

    # Run pitch sweep for all Rm
    theta_c_sim = np.zeros(len(mirror_list))
    err_theta_c_sim = np.zeros(len(mirror_list))
    for j in range(0,len(mirror_list)):
        print(f"Running sweep {j + 1:.0f}/{len(mirror_list):.0f}")
        [theta_c_sim[j],err_theta_c_sim[j]] = sim_mirror_sweep(mirror_list[j],v0,N,CONST) # find simulated critical pitch angle for current mirror/Rm [deg]

    # After sweep is completed, run single sim at critical pitch angle to check mu drift
    for k in range(0,len(mirror_list)):
        mode = "mu"
        state_check = sim_single_mirror(mirror_list[k],np.radians(mirror_list[k].theta_c + 5),v0,mode,CONST)
        mu_drift = check_mu_conservation(state_check,mirror_list[k],CONST)

        [epsilon,adiabatic] = invariance_check(v0,mirror_list[k],CONST) # adiabaticity condition for given mirror and v0

        # print simulation metrics
        print(f"Rm = {mirror_list[k].Rm:.1f} | theta_c analytic = {mirror_list[k].theta_c:.1f} deg | theta_c sim = {theta_c_sim[k]:.1f} deg | error = {err_theta_c_sim[k]:.2f}% | epsilon = {epsilon:.3f} | mu drift at theta_c = {mu_drift*100:.3f}%")

    # Plot simulated vs analytical critical pitch angles
    plot_crit_curve(Rm,theta_c_sim)


def sim_single_tokamak(tokamak,theta,v0,mode,cycles,CONST):
    """
    Run a single particle toroidal field simulation

    Parameters
    ----------
    tokamak : Tokamak object
        Magnetic field object belonging to the Tokamak class stored in tokamak.py

    theta : scalar
        Particle pitch angle in radians

    v0 : scalar
        Initial particle speed in meters per second

    mode : string
        "sweep" | "single" | "banana"
        Controls integration time and termination events
        sweep - t_end determined from estimated drift, no termination events
        single - termiantes after specified number of poloidal cycles through transit event
        banana - fixed t_end based on bounce period, no termination events

    cycles : scalar
        Number of poloidal cycles to integrate over, unitless

    CONST : named tuple, shape (2,)
        Named tuple containing the physical constants q (proton charge) in Coulombs and m (proton mass) in kilograms
        
    Returns
    -------
    state : OdeResult
        Solution object from solve_ivp
    """
        
    # Initial conditions
    v0_perp = v0 * np.sin(theta) # initial particle velocity perpendicular to B field [m/s]
    v0_par = v0 * np.cos(theta) # initial particle velocity parallel to B field [m/s]

    # Package inital state vector
    if mode == "banana":
        state0 = [tokamak.R0 + 0.3, 0.0, 0.0, # start slightly outboard of magnetic axis [m]
                  0.0, 0.0, v0] # all velocity perpendicular to B [m/s]
    else:
        state0 = [tokamak.R0, 0.0, 0.0, # initial position [m]
                  v0_perp, v0_par, 0.0] # initial velocity [m/s]
    
    # Integration time conditions computed based on cyclotron period and estimated drift speed
    omega_c = CONST.q * tokamak.B0 / CONST.m # cyclotron frequency [1/s]
    T_cyclotron = 2 * np.pi / omega_c # cyclotron period [s]
    max_step = T_cyclotron / 200 # 200 steps per gyration for max timestep [s]

    r_L = v0_perp / omega_c # Larmor radius [m]
    est_gradB = CONST.m * v0_perp**2 / (2 * CONST.q * tokamak.B0 * tokamak.R0) * (1 + r_L**2 / tokamak.R0**2) # estimated grad B drift [m/s]
    est_curvature = CONST.m * v0_par**2 / (CONST.q * tokamak.B0 * tokamak.R0) # estimated curvature drift [m/s]
    est_drift = est_gradB + est_curvature # total esimated drift [m/s]

    circ = 2 * np.pi * tokamak.R0 # tokamak circumference [m]
    T_toroidal = circ / max(abs(v0_par), 1e-6) # estimate toroidal period [s]

    if tokamak.q_safety is not None:
        T_poloidal = tokamak.q_safety * T_toroidal # estimate poloidal period [s]
    else: 
        T_poloidal = T_toroidal # if pure toroidal field there is no poloidal period, set as toroidal period

    # Define event to terminate integration after particle makes specific number of poloidal cycles
    def transit(t,state,tokamak,CONST):
        if t < 0.90 * cycles * T_poloidal: # don't allow event to trigger until 90% of predicted time passes
            return 1
        eps = 1e-10 # small offset to subtract from initial y position so that event doesn't trigger at start of integration
        return state[1] - (state0[1] - eps) # full transit is when y position passes intial y position
    
    transit.terminal = True # stop integration if particle completes a full transit
    transit.direction = 1 # stop integration when particle crosses initial x going positive

    if mode == "sweep":
        # adapative t_end runs simulation until drift covers 20 Larmor radii, with floor of 200 gyrations to ensure z drift dominates noise from individual gyrations
        t_end = max(20 * r_L / est_drift, 200 * T_cyclotron) # integration end time [s]
        t_span = [0,t_end] # integration time bounds [s]

        state = integrate.solve_ivp(accel, t_span, state0,
                                    method="RK45", max_step=max_step,
                                    rtol=1e-10, atol=1e-12,
                                    args=(tokamak,CONST))
    elif mode == "single":
        # 1.5x multiplier adds extra margin to make sure transit event is what ends integration
        t_end = 1.5* cycles * T_poloidal # integration end time [s]
        t_span = [0,t_end] # integration time bounds [s]

        state = integrate.solve_ivp(accel, t_span, state0,
                                    method="RK45", max_step=max_step,
                                    rtol=1e-10, atol=1e-12,
                                    args=(tokamak,CONST),
                                    events=[transit])
    elif mode == "banana":
        # Not currently implemented -- for use in milestone 3
        # Estimate bounce period from minor radius transit time
        # Run for 100 bounce periods to fill banana shape in RZ plot
        T_bounce_est = 4 * tokamak.r0 / v0 # estimate bounce period [s]
        t_end = 100 * T_bounce_est # integration end time [s]
        t_span = [0,t_end] # integration time bounds [s]

        state = integrate.solve_ivp(accel, t_span, state0,
                                    method="RK45", max_step=max_step,
                                    rtol=1e-10, atol=1e-12,
                                    args=(tokamak,CONST))
    
    return state


def sim_toroid_poloid_comp(toroid,poloid,v0,theta,CONST):
    """
    Run a single particle simulation for a pure toroidal field and combined field to evaluate drift suppression

    Parameters
    ----------
    toroid : Tokamak object
        Magnetic field object for the pure toroidal field

    poloid : Tokamak object
        Magnetic field object for the combined toroidal and poloidal field

    v0 : scalar
        Initial particle speed in meters per second

    theta : scalar
        Particle pitch angle in radians

    CONST : named tuple, shape (2,)
        Named tuple containing the physical constants q (proton charge) in Coulombs and m (proton mass) in kilograms
        
    Returns
    -------
    none
    """

    # Initial conditions
    v0_perp = v0 * np.sin(theta) # initial particle velocity perpendicular to B field [m/s]
    v0_par = v0 * np.cos(theta) # initial particle velocity parallel to B field [m/s]

    state0 = [poloid.R0, 0.0, 0.0, # initial position [m]
              v0_perp, v0_par, 0.0] # initial velocity [m/s]
    
    # Integration time conditions computed based on particle velocity and mirror length
    omega_c = CONST.q * poloid.B0 / CONST.m # cyclotron frequency [1/s]
    T_cyclotron = 2 * np.pi / omega_c # cyclotron period [s]
    max_step = T_cyclotron / 200 # 200 steps per gyration for max timestep [s]
    
    circ = 2 * np.pi * poloid.R0
    T_toroidal = circ / max(abs(v0_par), 1e-6)
        
    T_poloidal = poloid.q_safety * T_toroidal

    t_end = T_poloidal
    t_span = [0,t_end]


    state_tor = integrate.solve_ivp(accel, t_span, state0,
                                method="RK45", max_step=max_step,
                                rtol=1e-10, atol=1e-12,
                                args=(toroid,CONST))
    
    state_pol = integrate.solve_ivp(accel, t_span, state0,
                                method="RK45", max_step=max_step,
                                rtol=1e-10, atol=1e-12,
                                args=(poloid,CONST))
    
    plot_toroid_poloid_comp(state_tor,state_pol,poloid,theta)
    

def sim_tokamak_pitch_sweep(tokamak,v0,N,CONST):
    """
    Run simulations sweeping through pitch angles to compare to analytic drift predictions

    Parameters
    ----------
    tokamak : Tokamak object
        Magnetic field object belonging to the Tokamak class stored in tokamak.py

    v0 : scalar
        Initial particle speed in meters per second

    N : scalar
        Number of pitch angles to sweep through, unitless

    CONST : named tuple, shape (2,)
        Named tuple containing the physical constants q (proton charge) in Coulombs and m (proton mass) in kilograms
        
    Returns
    -------
    none
    """
        
    theta = np.radians(np.linspace(0,90,N)) # theta list to sweep

    # initialize empty arrays for results
    drift_sim = np.zeros(len(theta))
    drift_analytic = np.zeros(len(theta))
    err = np.zeros(len(theta))

    # run sweep
    for i in range(0,len(theta)):
        print(f"Running simulation {i+1}/{len(theta)}")
        mode = "sweep"
        cycles = 1
        state = sim_single_tokamak(tokamak,theta[i],v0,mode,cycles,CONST)
        drift_sim[i], drift_analytic[i], err[i] = check_drift_rate(state,tokamak,v0,theta[i],CONST)

    # plot results
    plot_drift(drift_sim,theta,err,tokamak,v0,CONST)


def sim_tokamak_q_sweep(q_list,R0,r0,B0,v0,theta,CONST):
    """
    Run tokamak simulations sweeping through a list of safety factors to validate linear relationship between q and z excursion

    Parameters
    ----------
    q_list : array-like, size (N,)
        List containing the safety factor values to sweep through

    R0 : scalar
        Major radius of tokamak in meters

    r0 : scalar
        Minor radius of tokamak in meters

    B0 : scalar
        Magnetic field strength at tokmaka major radius in Tesla

    v0 : scalar
        Initial particle speed in meters per second

    theta : scalar
        Particle pitch angle in radians

    CONST : named tuple, shape (2,)
        Named tuple containing the physical constants q (proton charge) in Coulombs and m (proton mass) in kilograms
        
    Returns
    -------
    none
    """
    
    z_ex = np.zeros(len(q_list)) # initialize empty array for excursion

    # run sweep
    for i in range(0,len(q_list)):
        tokamak = Tokamak(R0,r0,B0,q_list[i])
        mode = "single"
        cycles = 1
        state = sim_single_tokamak(tokamak,theta,v0,mode,cycles,CONST)
        z_ex[i] = np.max(state.y[2]) - np.min(state.y[2])
        
    # plot results
    plot_q_sweep(q_list,z_ex,tokamak,v0,theta,CONST)
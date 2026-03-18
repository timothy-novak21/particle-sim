import numpy as np
from collections import namedtuple

from fields.magnetic_mirror import Mirror
from fields.tokamak import Tokamak
from simulation import sim_single_mirror, sim_rm_sweep, sim_single_tokamak, sim_tokamak_pitch_sweep, sim_toroid_poloid_comp, sim_tokamak_q_sweep
from analysis.plotting import plot_single, plot_single_tokamak, plot_colormap


# Milestone 1
RUN_MIRROR_SWEEP = False
RUN_MIRROR_SINGLE = False

# Milestone 2
RUN_TOROIDAL_SWEEP = False
RUN_TOROIDAL_SINGLE = False
RUN_POLOIDAL_SINGLE = False
RUN_TOROIDAL_POLOIDAL_COMP = False
RUN_Q_SAFETY_SWEEP = False
RUN_TOKAMAK_COLOR = False

# Constants: q = proton charge [C], m = proton mass [kg]
PhysicalConstants = namedtuple("PhysicalConstants", ["q","m"])
CONST = PhysicalConstants(q=1.602e-19, m=1.673e-27)

# Define mirror ratio array and max B field
Bmax = 0.5 # B field strength at mirror [T]

# Define and print initial particle speed
KE_eV = 5e3 # physically relevant kinetic energy for a proton at fusion conditions [eV]
v0 = (CONST.q * KE_eV / (0.5 * CONST.m))**0.5 # solve for particle velocity based on desired KE [m/s]
# print(f"Initial particle kinetic energy: {KE_eV/1e3:.1f} keV, Initial particle velocity: {v0/1e6:.2f}e6 m/s")


if RUN_MIRROR_SWEEP:
    # Run full mirror pitch sweep for all Rm
    Rm = [1.5, 2, 3, 5, 10] # array of mirror ratios [-]
    N = 500 # number of pitch angles to sweep for mirror [-]
    sim_rm_sweep(Rm,Bmax,v0,N,CONST)

if RUN_MIRROR_SINGLE:
    # Set up and run single mirror sim
    length = 1 # half length of mirror [m]
    Rm_single = 2 # mirror ratio for single sim [-]
    mirror = Mirror(length,Rm_single,Bmax) # create mirror for single sim
    theta = np.radians(mirror.theta_c + 5) # run at theta slightly higher than critical pitch [rad]
    mode = "single"
    state = sim_single_mirror(mirror,theta,v0,mode,CONST)
    plot_single(state,mirror,theta) # plot single


# Set up toroidal field
R0 = 3 # toroidal field major radius [m]
r0 = 0.9 # toroidal field minor radius [m]
B0 = Bmax # field strength along toroid centerline [T]
toroid = Tokamak(R0,r0,B0) # initialize tokamak object with pure toroidal field
q_safety = 2.0 # safety factor for poloidal field strength [-]
poloid = Tokamak(R0,r0,B0,q_safety) # initialize tokamak object with poloidal correction

if RUN_TOROIDAL_SWEEP:
    # Run toroidal pitch sweep
    N = 101 # number of pitch angles to sweep for toroidal [-]
    sim_tokamak_pitch_sweep(toroid,v0,N,CONST)

if RUN_TOROIDAL_SINGLE:
    # Set up and run single sim with pure toroidal field
    theta = np.radians(20) # pitch angle [rad]
    mode = "single"
    cycles = 1
    state = sim_single_tokamak(toroid,theta,v0,mode,cycles,CONST)
    plot_single_tokamak(state,toroid,theta)

if RUN_POLOIDAL_SINGLE:
    # Set up and run single sim with poloidal correction
    theta = np.radians(20) # pitch angle [rad]
    mode = "single"
    cycles = 1
    state = sim_single_tokamak(poloid,theta,v0,mode,cycles,CONST)
    plot_single_tokamak(state,poloid,theta)

if RUN_TOROIDAL_POLOIDAL_COMP:
    theta = np.radians(20) # pitch angle [rad]
    sim_toroid_poloid_comp(toroid,poloid,v0,theta,CONST)

if RUN_Q_SAFETY_SWEEP:
    q_list = [1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10] # list of safety factors [-]
    theta = np.radians(20) # pitch angle [rad]
    sim_tokamak_q_sweep(q_list,R0,r0,B0,v0,theta,CONST)

if RUN_TOKAMAK_COLOR:
    # Set up and run tokamak sim with 3d color plot
    theta = np.radians(20) # pitch angle [rad]
    mode = "single"
    cycles = 1
    state = sim_single_tokamak(poloid,theta,v0,mode,cycles,CONST)
    plot_colormap(state,poloid,theta)

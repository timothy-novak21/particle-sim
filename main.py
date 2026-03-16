import numpy as np
from collections import namedtuple

from fields.magnetic_mirror import Mirror
from simulation import sim_single, sim_rm_sweep
from analysis.plotting import plot_single


# Constants: q = proton charge [C], m = proton mass [kg]
PhysicalConstants = namedtuple("PhysicalConstants", ["q","m"])
CONST = PhysicalConstants(q=1.602e-19, m=1.673e-27)

# Define mirror ratio array and max B field
Rm = [1.5, 2, 3, 5, 10] # array of mirror ratios [-]
Bmax = 0.5 # B field strength at mirror [T]

# Define and print initial particle speed
KE_eV = 5e3 # physically relevant kinetic energy for a proton at fusion conditions [eV]
v0 = (CONST.q * KE_eV / (0.5 * CONST.m))**0.5 # solve for particle velocity based on desired KE [m/s]
print(f"Initial particle kinetic energy: {KE_eV/1e3:.1f} keV, Initial particle velocity: {v0/1e6:.2f}e6 m/s")

N = 500 # number of pitch angles to sweep [-]

RUN_SWEEP = True
RUN_SINGLE = False

if RUN_SWEEP:
    # Run full pitch sweep for all Rm
    sim_rm_sweep(Rm,Bmax,v0,N,CONST)

if RUN_SINGLE:
    # Set up and run single sim
    length = 1 # half length of mirror [m]
    Rm_single = 2 # mirror ratio for single sim [-]
    mirror = Mirror(length,Rm_single,Bmax) # create mirror for single sim
    theta = np.radians(mirror.theta_c + 5) # run at theta slightly higher than critical pitch [rad]
    state = sim_single(mirror,theta,v0,"single",CONST)
    plot_single(state,mirror,theta) # plot single

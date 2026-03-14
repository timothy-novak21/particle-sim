import numpy as np

def invariance_check(v0,mirror,CONST):
    """
    Calculate the adiabaticity parameter and check invariance for a magnetic mirror

    Parameters
    ----------
    v0 : scalar
        Initial particle speed in meters per second

    mirror : Mirror object
        Magnetic mirror object belonging to the Mirror class stored in B_field.py

    CONST : named tuple, shape (2,)
        Named tuple containing the physical constants q (proton charge) in Coulombs and m (proton mass) in kilograms
        
    Returns
    -------
    epsilon : scalar
        Adiabaticity parameter that is used to determine if the B field is adiabatically invariant, unitless

    invariant : boolean
        Stores True if the B field is invariant (epsilon << 1) and False if the field is NOT invariant
    """

    invariant = False

    r_L = (CONST.m * v0) / (CONST.q * mirror.Bmax) # Larmor radius at mirror point [m]
    dBz_dz = 2 * (mirror.Bmax - mirror.Bmin) / mirror.length # field gradient magnitude at mirror boundary
    epsilon = r_L * dBz_dz / mirror.Bmax # adiabaticity parameter [-]

    # eps << 1 for adiabaticity to hold
    if epsilon < 0.01:
        invariant = True

    return [epsilon,invariant]
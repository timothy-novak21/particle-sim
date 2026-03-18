import numpy as np


class Tokamak:
    """
    Class to store variables and functions for a toroidal magnetic field

    Variables
    ---------
    R0 : scalar
        Major radius of torus in meters

    r0 : scalar
        Minor radius of torus in meters
    
    B0 : scalar
        Magnetic field strength of toroidal field at R0 in Tesla

    q_safety : scalar or None
        Safety factor defining poloidal field strength, unitless
        If None, tokamak object is a pure toroidal field with no poloidal component
    """

    def __init__(self,R0,r0,B0,q_safety=None):
        self.R0 = R0 # toroidal field major radius [m]
        self.r0 = r0 # toroidal field minor radius [m]
        self.B0 = B0 # B field strength at R0 [T]
        self.q_safety = q_safety # safety factor [-]


    def vector(self,pos):
        """
        Compute the magnetic field vector at a given position in a tokamak field

        Parameters
        ----------
        pos : array-like, shape (3,)
            Cartesian position [x,y,z] in meters.
        
        Returns
        -------
        Bvec : np.array, shape(3,)
            Magnetic field vector [Bx,By,Bz] in Tesla
        """

        x,y,z = pos # unpack cartesian coordinates [m]
        eps = 1e-10 # epsilon to prevent division by zero during cartesian conversion
        R = (x**2 + y**2)**0.5 + eps # radial position [m]

        B_tor = self.B0 * self.R0 / R # magnitude of toroidal field at given radial position [T]

        # Project toroidal field into cartesian coordinates
        Bx = B_tor * (-y / R)
        By = B_tor * (x / R)
        Bz = 0

        # Add poloidal field (if q_safety exists)
        if self.q_safety is not None:
            r_minor = ((R - self.R0)**2 + z**2)**0.5 + eps # radial offset inside torus [m]
            
            B_pol = (r_minor * B_tor) / (self.R0 * self.q_safety) # magnitude of poloidal field at given position [T]

            # Add poloidal field correction
            Bx += -B_pol * (z / r_minor) * (x / R)
            By += -B_pol * (z / r_minor) * (y / R)
            Bz += B_pol * ((R - self.R0) / r_minor)

        Bvec = np.array([Bx, # B-field x-component [T]
                         By, # B-field y-component [T]
                         Bz]) # B-field z-component [T]
        
        return Bvec
    

    def field_check(self):
        """
        Check that B field equations are set up properly by calculated B field strength at various locations in toroidal field

        For use as a manual validation tool during development
        Call toroid.field_check() to verify field equations

        note: currently field_check() is only set up to check pure toroidal fields
        todo: update to allow validation of combined fields

        Parameters
        ----------
        none
        
        Returns
        -------
        none
        """

        # All position sets in cartesian [m]
        # Positive x, central z position set
        pos_inner_posx_z0 = [self.R0 - self.r0, 0.0, 0.0] # inner torus wall
        pos_outer_posx_z0 = [self.R0 + self.r0, 0.0, 0.0] # outer torus wall
        pos_central_posx_z0 = [self.R0, 0.0, 0.0] # center of torus

        # Positive y, central z position set
        pos_inner_posy_z0 = [0.0, self.R0 - self.r0, 0.0]
        pos_outer_posy_z0 = [0.0, self.R0 + self.r0, 0.0]
        pos_central_posy_z0 = [0.0, self.R0, 0.0]

        # Negative x, central z position set
        pos_inner_negx_z0 = [-(self.R0 - self.r0), 0.0, 0.0]
        pos_outer_negx_z0 = [-(self.R0 + self.r0), 0.0, 0.0]
        pos_central_negx_z0 = [-self.R0, 0.0, 0.0]

        # Negative y, central z position set
        pos_inner_negy_z0 = [0.0, -(self.R0 - self.r0), 0.0]
        pos_outer_negy_z0 = [0.0, -(self.R0 + self.r0), 0.0]
        pos_central_negy_z0 = [0.0, -self.R0, 0.0]

        # Positive x, shifting z position set
        pos_central_posx_posz = [self.R0, 0.0, self.r0] # postive z / upper torus wall
        pos_central_posx_negz = [self.R0, 0.0, -self.r0] # negative z / lower torus wall

        # Package positions
        pos = [pos_inner_posx_z0,pos_outer_posx_z0,pos_central_posx_z0,
               pos_inner_posy_z0,pos_outer_posy_z0,pos_central_posy_z0,
               pos_inner_negx_z0,pos_outer_negx_z0,pos_central_negx_z0,
               pos_inner_negy_z0,pos_outer_negy_z0,pos_central_negy_z0,
               pos_central_posx_posz,pos_central_posx_negz]

        # Intiailize empty Bmag and exp_B arrays
        Bmag = np.zeros(len(pos))
        exp_B = np.zeros(len(pos))
        for i in range(0,len(pos)):
            Bvec = self.vector(pos[i]) # find B field vector at given position [T]
            Bmag[i] = np.linalg.norm(Bvec) # find magnitude of B field at given position [T]
            exp_B[i] = self.B0 * self.R0 / (pos[i][0]**2 + pos[i][1]**2)**0.5 # calc expected B field strength

        # Print calculated vs expected B field
        print(f"B at inner wall, positive x, z=0: {Bmag[0]:.4f} T (expected {exp_B[0]:.4f})")
        print(f"B at outer wall, positive x, z=0: {Bmag[1]:.4f} T (expected {exp_B[1]:.4f})")
        print(f"B at center, positive x, z=0: {Bmag[2]:.4f} T (expected {exp_B[2]:.4f})")
        print(f"-----------------------------------------------------------")

        print(f"B at inner wall, positive y, z=0: {Bmag[3]:.4f} T (expected {exp_B[3]:.4f})")
        print(f"B at outer wall, positive y, z=0: {Bmag[4]:.4f} T (expected {exp_B[4]:.4f})")
        print(f"B at center, positive y, z=0: {Bmag[5]:.4f} T (expected {exp_B[5]:.4f})")
        print(f"-----------------------------------------------------------")

        print(f"B at inner wall, negative x, z=0: {Bmag[6]:.4f} T (expected {exp_B[6]:.4f})")
        print(f"B at outer wall, negative x, z=0: {Bmag[7]:.4f} T (expected {exp_B[7]:.4f})")
        print(f"B at center, negative x, z=0: {Bmag[8]:.4f} T (expected {exp_B[8]:.4f})")
        print(f"-----------------------------------------------------------")

        print(f"B at inner wall, negative y, z=0: {Bmag[9]:.4f} T (expected {exp_B[9]:.4f})")
        print(f"B at outer wall, negative y, z=0: {Bmag[10]:.4f} T (expected {exp_B[10]:.4f})")
        print(f"B at center, negative y, z=0: {Bmag[11]:.4f} T (expected {exp_B[11]:.4f})")
        print(f"-----------------------------------------------------------")

        print(f"B at center, positive x, z=0: {Bmag[2]:.4f} T (expected {exp_B[2]:.4f})")
        print(f"B at center, positive x, positive z: {Bmag[12]:.4f} T (expected {exp_B[12]:.4f})")
        print(f"B at center, positive x, positive z: {Bmag[13]:.4f} T (expected {exp_B[13]:.4f})")
        print(f"-----------------------------------------------------------")

        
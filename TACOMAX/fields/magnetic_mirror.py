import numpy as np


class Mirror:
    """
    Class to store variables and functions for a magnetic mirror

    Variables
    ---------
    length : scalar
        Half length of magnetic mirror in meters

    Rm : scalar
        Mirror ratio, equal to Bmax divided by Bmin, unitless
    
    Bmax : scalar
        Maximum magnetic field strength, located at mirror boundary in Tesla

    Bmin : scalar
        Minimum magnetic field strength, located at mirror center in Tesla

    theta_c : scalar
        Critical pitch angle in degrees, all angles less than this will escape the mirror
    """


    def __init__(self,length,Rm,Bmax):
        self.length = length # Mirror half length [m]
        self.Rm = Rm # Mirror ratio [-]
        self.Bmax = Bmax # Magnetic field strength at mirror boundary [T]
        self.Bmin = Bmax / Rm # Magnetic field strength at mirror center [T]
        self.theta_c = np.degrees(np.arcsin(1/self.Rm**0.5)) # critical pitch angle [deg]


    def vector(self,pos):
        """
        Compute the magnetic field vector at a given position for a magnetic mirror

        Bz can be found via Biot-Savart law, however a simpler approximation can be 
        found in the form Bz = A + C*z^2, where A = Bmin and C = (Bmax - Bmin)/L^2

        Maxwell's equations require div(B) = 0. In cylindrical coordinates this becomes
        (1/r)(d(r*Br)/dr) + dBz/dz = 0, which integrates to Br = (-r/2)*(dBz/dz)

        For the given Bz, this derivative becomes dBz/dz = 2*(Bmax - Bmin)*z/L^2

        Parameters
        ----------
        pos : array-like, shape (3,)
            Cartesian position [x,y,z] in meters.
        
        Returns
        -------
        Bvec : np.array, shape(3,)
            Magnetic field vector [Bx,By,Bz] in Tesla
        """

        x,y,z = pos # unpack cartesian coordinates from position [m]
        eps = 1e-10 # epsilon to prevent division by zero during cartesian conversion
        r = (x**2 + y**2)**0.5 + eps # radial position [m]

        # Analytic mirror field approximation
        Bz = self.Bmin + (self.Bmax - self.Bmin) * (z / self.length)**2

        # Radial component from div(B) = 0
        dBz_dz = 2 * (self.Bmax - self.Bmin) * z / self.length**2
        Br = (-r/2) * dBz_dz

        # Project radial component into cartesian x and y
        Bx = Br * (x / r)
        By = Br * (y / r)

        # Combine magnetic field components into vector
        Bvec = np.array([Bx, # B-field x-component [T]
                         By, # B-field y-component [T]
                         Bz]) # B-field z-component [T]

        return Bvec
    

    def field_check(self):
        """
        Check that B field equations are set up properly by calculated B field strength at mirror center and boundary

        For use as a manual validation tool during development
        Call mirror.field_check() to verify field equations

        Parameters
        ----------
        none
        
        Returns
        -------
        none
        """
        
        pos_center = np.array([0.0, 0.0, 0.0]) # cartesian position at center of mirror [m]
        pos_boundary = np.array([0.0, 0.0, self.length]) # cartesian position at boundary of mirror [m]

        Bvec_center = self.vector(pos_center) # B field vector at mirror center [T]
        Bvec_boundary = self.vector(pos_boundary) # B field vector at mirror boundary [T]

        Bmag_center = np.linalg.norm(Bvec_center) # B field magnitude at mirror center [T]
        Bmag_boundary = np.linalg.norm(Bvec_boundary) # B field magnitude at mirror boundary [T]

        # Display calculated vs expected B field values
        print(f"B at center: {Bmag_center:.4f} T (expected {self.Bmin:.4f})")
        print(f"B at mirror boundary: {Bmag_boundary:.4f} T (expected {self.Bmax:.4f})")
        print(f"Actual mirror ratio: {Bmag_boundary/Bmag_center:.4f} (expected {self.Rm:.4f})")
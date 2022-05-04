import numpy as np
import matplotlib.pyplot as plt

# Constants
K2 = 0.425895581    # a derived constant [s^3/m^3]
rho_sea = 1.225     # air density at sea-level [kg/m^3]
rho_cr = 1.17315    # air density at cruise altitude [kg/m^3]
h_cr = 450          # cruise altitude [m]

def rotor_solidity(B: int, c: float, R: float):
    """
    This function calculates the rotor solidity

    Parameters
    ----------
    B: int
        Number of blades out of which the rotor consists [m]

    c: float
        Average blade chord [m]

    R: float
        Blade radius [m]

    Returns
    -------
    s: float
        rotor solidity [-]
    """

    s = (B*c)/(np.pi*R)
    return s

def tip_velocity():
    return 0

def vortex_noise(B: int, c: float, R: float, rho: float, Delta_S: float, N: int, T: float, disk_loading: float):
    """
    This function calculates the vortex noise produced by the vehicle [in dB]

    Parameters
    ----------
    B: int
        Number of blades out of which the rotor consists [m]

    c: float
        Average blade chord [m]

    R: float
        Blade radius [m]

    rho: float
        Density at altitude [kg/m^3]

    Delta_S: float
        Distance to observer [m]

    N: float
        Number of rotors present on vehicle [-]

    T: float
        Thrust required per engine [N]

    disk_loading: float
        The thrust per rotor disk area [N/m^2]

    Returns
    -------
    SPL_vortex: float
        vortex noise [dB]

    """
    s = rotor_solidity(B, c, R)
    V_T = tip_velocity()
    SPL_vortex = 20 * np.log10(K2*(V_T/(rho*Delta_S)) * np.sqrt((N*T/s) * disk_loading))
    return SPL_vortex
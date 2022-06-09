import numpy as np

# Constants
g = 9.80665         # gravitation acceleration [m/s^2]
K2 = 0.425895581    # a derived constant [s^3/m^3]
rho_cr = 1.225     # air density at sea-level [kg/m^3]
# rho_cr = 1.17315    # air density at cruise altitude [kg/m^3]
h_cr = 450       # cruise altitude [m]
Cl = 1.0
disk_loading = 500
c = 0.1
B = 3

# Functions for the noise calculations
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
        Rotor solidity [-]
    """

    s = (B*c)/(np.pi*R)
    return s

def tip_velocity(B: int, c: float, R: float, rho: float, Thrust: float, Cl: float):
    """
    This function calculates the tip velocity of the rotor's propeller

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

    Thrust: float
        Thrust provided/required by ONE engine [N]

    Cl: float
        Blade mean lift coefficient [-]

    Returns
    -------
    V_T: float
        Propeller tip velocity [m/s]

    """
    V_T = np.sqrt((6*Thrust)/(rho * B * c * R * Cl))
    return V_T

def vortex_noise(N: int, T_total: float, R: float):
    B = 3
    c = 0.1
    rho = rho_cr
    Delta_s = h_cr
    Cl = 1.0
    # disk_loading = 500
    Thrust_engine = T_total / N
    # print(Thrust_engine)
    A_rotor = np.pi*R** 2
    disk_loading = Thrust_engine / A_rotor
    # R = np.sqrt(A_rotor / np.pi)
    s = rotor_solidity(B, c, R)
    V_t = np.sqrt((6*np.sqrt(disk_loading*np.pi*Thrust_engine)/(rho*B*c*Cl)))
    # print(s, V_t, A_rotor, R, Thrust_engine)
    SPL_vortex = K2*(V_t/(rho*Delta_s)) * np.sqrt((N*Thrust_engine/s) * disk_loading)
    return SPL_vortex

def SPL_calc(noise_intensity):
    SPL = 20*np.log10(noise_intensity)
    return SPL

intensity_2V = vortex_noise(4, 2955*4, 1.56376)
intensity_2H = vortex_noise(1, 2654.25, 1.3)
total_intensity_2 = intensity_2H + intensity_2V
SPL_2 = SPL_calc(total_intensity_2)
print(SPL_2)

intensity_4V = vortex_noise(4, 1623.78*g, 1.109)
intensity_4H = vortex_noise(2, 1401.324556721424, 0.944516)
total_intensity_4 = intensity_4H + intensity_4V
SPL_4 = SPL_calc(total_intensity_4)
print(SPL_4)
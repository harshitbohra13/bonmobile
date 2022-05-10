import numpy as np
import matplotlib.pyplot as plt

# Constants
g = 9.80665         # gravitational acceleration
K2 = 0.425895581    # a derived constant [s^3/m^3]
rho_sea = 1.225     # air density at sea-level [kg/m^3]
rho_cr = 1.17315    # air density at cruise altitude [kg/m^3]
h_cr = 450          # cruise altitude [m]
d_to_build = 4.6    # closest distance to a building the vehicle is allowed to get [m]

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

def vortex_noise(B: int, c: float, R: float, rho: float, Delta_S: float, N: int, Thrust: float, Cl: float, disk_loading: float):
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

    Thrust: float
        Thrust required per engine [N]

    Cl: float
        Blade mean lift coefficient [-]

    disk_loading: float
        The thrust per rotor disk area [N/m^2]

    Returns
    -------
    SPL_vortex: float
        vortex noise [dB]

    """
    s = rotor_solidity(B, c, R)
    V_T = tip_velocity(B, c, R, rho, Thrust, Cl)
    SPL_vortex = 20 * np.log10(K2*(V_T/(rho*Delta_S)) * np.sqrt((N*Thrust/s) * disk_loading))
    return SPL_vortex

# Providing arrays with values for each concept
mass = np.array([1452.483, 1599.3759, 1520.29, 1623.78, 1550.758])
Weight = mass * g
B = np.array([3, 3, 3, 3, 3])
c = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
N = np.array([4, 4, 8, 8, 12])
Thrust_engine = Weight/N
Cl = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
disk_loading = np.array([500, 500, 500, 500, 500])
A = []
R = []
for i in range(len(disk_loading)):
    A.append(Thrust_engine[i]/disk_loading[i])
    R.append(np.sqrt(A[i]/np.pi))

SPL_vortex = []
s = []
for i in range(5):
    SPL_vortex.append(vortex_noise(B[i], c[i], R[i], rho_cr, h_cr, N[i], Thrust_engine[i], Cl[i], disk_loading[i]))
    s.append(rotor_solidity(B[i], c[i], R[i]))

# print(disk_loading)

sens_disk_loading_list = [*range(100, 500, 50)]

SPL_vortex_list_rho1 = []
SPL_vortex_list_rho2 = []
for i in range(len(sens_disk_loading_list)):
    SPL_vortex_list_rho1.append(vortex_noise(4, 0.1, 0.5, rho_cr, h_cr, 4, 300, 0.8, sens_disk_loading_list[i]))
    SPL_vortex_list_rho2.append(vortex_noise(4, 0.1, 0.5, rho_sea, h_cr/2, 4, 300, 0.8, sens_disk_loading_list[i]))

plt.figure(1)
plt.title("Noise sensitivity to disk-loading")
plt.plot(sens_disk_loading_list,SPL_vortex_list_rho1, label="cruise altitude")
plt.plot(sens_disk_loading_list,SPL_vortex_list_rho2, label="sea-level")
plt.ylim([30, 60])
plt.grid()
plt.legend()


sens_NBlades_list = [*range(1, 9, 1)]

SPL_vortex_list_B = []
for i in range(len(sens_NBlades_list)):
    SPL_vortex_list_B.append(vortex_noise(sens_NBlades_list[i], 0.1, 0.5, rho_cr, h_cr, 4, 300, 0.8, 200))

plt.figure(2)
plt.title("Noise sensitivity to amount of rotor blades")
plt.plot(sens_NBlades_list, SPL_vortex_list_B)
# plt.yscale("log")
plt.ylim([30, 50])
plt.grid()


sens_Cl_list = np.arange(0.1, 1.6, 0.1)

SPL_vortex_list_Cl = []
for i in range(len(sens_Cl_list)):
    SPL_vortex_list_Cl.append(vortex_noise(4, 0.1, 0.5, rho_cr, h_cr, 4, 300, sens_Cl_list[i], 200))

plt.figure(3)
plt.title("Noise sensitivity to mean blade lift coefficient")
plt.plot(sens_Cl_list, SPL_vortex_list_Cl)
# plt.yscale("log")
plt.ylim([30, 50])
plt.grid()
plt.show()

sens_N_list = np.arange(1, 15, 1)

SPL_vortex_list_N = []
for i in range(len(sens_N_list)):
    SPL_vortex_list_N.append(vortex_noise(4, 0.1, 0.5, rho_cr, h_cr, sens_N_list[i], 1200/sens_N_list[i], 0.8, 200))

plt.figure(4)
plt.title("Noise sensitivity to amount of rotors")
plt.plot(sens_N_list, SPL_vortex_list_N)
# plt.yscale("log")
plt.ylim([30, 50])
plt.grid()
plt.show()

print("Noise produced by each vehicle:")
print("Noise vehicle 1 = " + str(round(SPL_vortex[0], 2)) + " dB")
print("Noise vehicle 2 = " + str(round(SPL_vortex[1], 2)) + " dB")
print("Noise vehicle 3 = " + str(round(SPL_vortex[2], 2)) + " dB")
print("Noise vehicle 4 = " + str(round(SPL_vortex[3], 2)) + " dB")
print("Noise vehicle 5 = " + str(round(SPL_vortex[4], 2)) + " dB")
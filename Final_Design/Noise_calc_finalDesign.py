import numpy as np
import matplotlib.pyplot as plt

# Constants
g = 9.80665         # gravitational acceleration
K2 = 0.425895581    # a derived constant [s^3/m^3]
rho_sea = 1.225     # air density at sea-level [kg/m^3]
rho_cr = 1.17315    # air density at cruise altitude [kg/m^3]
h_cr = 450          # cruise altitude [m]
d_req = 30.48       # distance set by requirement [m]

def addlabels(x,y,error):
    for i in range(len(x)):
        plt.text(i, y[i]/2, str(round(y[i],0)) + " +- " +  str(round(error[i], 1)) , ha = 'center')

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
mass = 1450
Weight = mass * g
B = 2
c = 0.092
N = 12
Thrust_engine = Weight/N
Cl = 1.0
R = 0.8
A = np.pi * R**2
disk_loading = Thrust_engine/A

SPL_vortex_cruise = vortex_noise(B, c, R, rho_cr, h_cr, N, Thrust_engine, Cl, disk_loading)
SPL_vortex_req = vortex_noise(B, c, R, rho_cr, d_req, N, Thrust_engine, Cl, disk_loading)
s = rotor_solidity(B, c, R)

# sens_disk_loading_list = [*range(1, 2000, 50)]
#
# SPL_vortex_list_rho1 = []
# SPL_vortex_list_rho2 = []
# for i in range(len(sens_disk_loading_list)):
#     SPL_vortex_list_rho1.append(vortex_noise(4, 0.1, 0.5, rho_cr, h_cr, 4, 300, 0.8, sens_disk_loading_list[i]))
#     SPL_vortex_list_rho2.append(vortex_noise(4, 0.1, 0.5, rho_sea, h_cr, 4, 300, 0.8, sens_disk_loading_list[i]))
#
# plt.figure(1)
# plt.title("Noise sensitivity to disk-loading")
# plt.plot(sens_disk_loading_list,SPL_vortex_list_rho1, label="cruise altitude")
# plt.plot(sens_disk_loading_list,SPL_vortex_list_rho2, label="sea-level")
# plt.ylim([20, 50])
# plt.grid()
# plt.legend()
#
#
# sens_NBlades_list = [*range(1, 14, 1)]
#
# SPL_vortex_list_B = []
# for i in range(len(sens_NBlades_list)):
#     SPL_vortex_list_B.append(vortex_noise(sens_NBlades_list[i], 0.1, 0.5, rho_cr, h_cr, 4, 300, 0.8, 200))
#
# plt.figure(2)
# plt.title("Noise sensitivity to amount of rotor blades")
# plt.plot(sens_NBlades_list, SPL_vortex_list_B)
# plt.ylim([10, 50])
# plt.grid()
#
#
# sens_Cl_list = np.arange(0.1, 1.6, 0.1)
#
# SPL_vortex_list_Cl = []
# for i in range(len(sens_Cl_list)):
#     SPL_vortex_list_Cl.append(vortex_noise(4, 0.1, 0.5, rho_cr, h_cr, 4, 300, sens_Cl_list[i], 200))
#
# plt.figure(3)
# plt.title("Noise sensitivity to mean blade lift coefficient")
# plt.plot(sens_Cl_list, SPL_vortex_list_Cl)
# plt.ylim([30, 50])
# plt.grid()
# plt.show()
#
# sens_N_list = np.arange(1, 15, 1)
#
# SPL_vortex_list_N = []
# for i in range(len(sens_N_list)):
#     SPL_vortex_list_N.append(vortex_noise(4, 0.1, 0.5, rho_cr, h_cr, sens_N_list[i], 1200/sens_N_list[i], 0.8, 200))
#
# plt.figure(4)
# plt.title("Noise sensitivity to amount of rotors")
# plt.plot(sens_N_list, SPL_vortex_list_N)
# # plt.yscale("log")
# plt.ylim([30, 50])
# plt.grid()


print("Noise produced by the vehicle at cruise altitude =", str(round(SPL_vortex_cruise, 2)), "dB")
print("Noise produced by the vehicle at 100 ft =", str(round(SPL_vortex_req, 2)), "dB")


mass_list = np.arange(1, 2500, 50)
thrust_list = [mass_list[i]*g/N for i in range(len(mass_list))]
noise_list_30 = [vortex_noise(B, c, R, rho_cr, d_req, N, thrust_list[i], Cl, thrust_list[i]/A) for i in range(len(thrust_list))]
noise_list_450 = [vortex_noise(B, c, R, rho_cr, h_cr, N, thrust_list[i], Cl, thrust_list[i]/A) for i in range(len(thrust_list))]

plt.figure(5)
plt.fill_between(mass_list, noise_list_30, noise_list_450, alpha=0.2)
plt.plot(mass_list, noise_list_30, c='blue', linestyle='--', label=f"$\Delta S$ = 30.48 m")
plt.plot(mass_list, noise_list_450, c='orange', linestyle='--', label=f"$\Delta S$ = 450 m")
plt.plot()
plt.xlabel(f"$m$ [kg]")
plt.ylabel(f"$SPL$ [dB]")
plt.legend()


distance_list = np.arange(0.1, 700, 10)
noise_list_plus10mass = [vortex_noise(B, c, R, rho_cr, distance_list[i], N, mass*1.1*g/N, Cl, (mass*1.1*g/N)/A) for i in range(len(distance_list))]
noise_list_minus10mass = [vortex_noise(B, c, R, rho_cr, distance_list[i], N, mass*0.9*g/N, Cl, (mass*0.9*g/N)/A) for i in range(len(distance_list))]

plt.figure(6)
plt.fill_between(distance_list, noise_list_plus10mass, noise_list_minus10mass, alpha=0.2)

plt.plot(distance_list, noise_list_plus10mass, c='blue', linestyle='--', label=f"+10% margin mass")
plt.plot(distance_list, noise_list_minus10mass, c='orange', linestyle='--', label=f"-10% margin mass")
# plt.axvline(x=450, color='r', linestyle='--', label="Cruise altitude")
# plt.axvline(x=30.48, color='g', linestyle='--', label="Requirement distance (100 ft)")
plt.axhline(y=70, color='r', linestyle='--', linewidth=0.7, label="Sound requirement 100 ft")
plt.plot(d_req, SPL_vortex_req, "o", label="Vehicle noise (100 ft)")
plt.plot(h_cr, SPL_vortex_cruise, "o", label="Vehicle noise (cruise altitude)")
plt.xlabel(f"$\Delta S$ [m]")
plt.ylabel(f"$SPL$ [dB]")
plt.xlim(0, 550)
plt.ylim(60, 90)
plt.legend()


# Define labels, positions, bar heights and error bar heights
labels = [f'$\Delta S = 30.48 m$', f'$\Delta S = 450 m$']
error_req = vortex_noise(B, c, R, rho_sea, d_req, N, mass*1.1*g/N, Cl, (mass*1.1*g/N)/A) - vortex_noise(B, c, R, rho_cr, d_req, N, mass*1.0*g/N, Cl, (mass*g/N)/A)
error_cr = vortex_noise(B, c, R, rho_cr, h_cr, N, mass*1.1*g/N, Cl, (mass*1.1*g/N)/A) - vortex_noise(B, c, R, rho_cr, h_cr, N, mass*1.0*g/N, Cl, (mass*g/N)/A)
x_pos = np.arange(len(labels))
noise_val = [SPL_vortex_req, SPL_vortex_cruise]
error = [error_req, error_cr]
fig, ax = plt.subplots()
ax.bar(x_pos, noise_val,
       yerr=error,
       align='center',
       alpha=0.6,
       ecolor='black',
       capsize=15)
addlabels(x_pos, noise_val, error)
ax.set_ylabel(f'$SPL$ [dB]')
ax.set_xticks(x_pos)
ax.set_xticklabels(labels)
ax.yaxis.grid(True)


plt.show()
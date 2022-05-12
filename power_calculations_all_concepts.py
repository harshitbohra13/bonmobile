import numpy as np
import drag_calc as drag

# Constants
g = 9.80665                                 # gravitation acceleration [m/s^2]
rho = 1.225                                 # density at sea-level [kg/m^3]
kmh_to_ms = 3.6                             # conversion factor from km/h to m/s [-]

# Aircraft data
mass_f = 500                                # fuselage mass [kg]
W_to = 1490 * g                             # Take-off weight CONCEPT 2 [N]
W_D = W_to                                  #
FOM = 0.75                                  # figure of merit [-]
CD0 = 0.3                                   # zero-lift drag [-]
S = 7.0                                     # surface area fuselage [m^2]
alpha = 5                                   #
C_D = 0.008                                 #
C_L = 0.5                                   # Lift coefficient [-]
S_wing = 15                                 # Surface area [-]
n_ult = 3.75                                # Ultimate load factor
AR = 8.3                                    # Aspect ratio [-]
b = np.sqrt(AR*S_wing)                      # Wing span [m]
t_c = 0.12                                  # Thickness-to-chord-ratio [-]
taper_ratio = 0.7                           # Wing taper ratio [-]
f_lambda = (1+ 2* taper_ratio)/(1+taper_ratio)
sweep = 0                                   # Wing sweep [deg]
specific_material_weight = 2700*g           # Specific weight of the used material [N/m^3]
t_ss = 0.004 * (1 + np.sqrt(W_to/(10**6)))
w_sigma = 40 * (1+1.1*(W_to/(10**6))**(-1/4))*10**(-6)
Weight_WING = 0.06 * w_sigma * n_ult * W_D * b**3 / S * f_lambda / \
              (t_c * np.cos(np.radians(sweep))**2) + specific_material_weight*t_ss*S


# Velocity data/requirements
V_descend = 7.0                             # descend velocity [m/s]
V_climb = 7.0                               # climb velocity [m/s]
V_cruise = 100/kmh_to_ms                    # cruise velocity [m/s]

# Time data
t_climb = (450 - 30.5)/V_climb              # time to climb [s]
t_descend = t_climb                         # time to descend [s]
t_cruise = 20/100 * 3600 + 120              # cruise time +120 seconds for acceleration and 20 km in one direction [s]
t_hover = 60 * 2                            # hovering time [s] -> hovering appears twice

battery_efficiency = 0.85                   # battery efficiency [-]
battery_density = 250                       # battery density [Wh/kg]

# ---concept 1---
# rotors data
blades_number = 3 * np.ones(5)              # number of blades per rotor [-]
rotors_number = np.array([4, 4, 8, 8, 12])  # number of rotors per concept [-]
disk_loading = 500 * np.ones(5)             # pre-defined constant disk-loading [N/m^2]

# Circular beam
structure_penalty = np.array([1, 1, 1.2, 1, 1.5])       # structure penalty for additional mass
structure_length = np.array([2, 1, 2, 4, 2])            # length of the beam [m]
structure_radius = 0.075                                # outer radius of the beam [m]
structure_inner_radius = 0.070                          # outer radius of the beam
structure_Area = np.pi * structure_radius ** 2 - \
                 np.pi * structure_inner_radius ** 2    # area of the beam
structure_density = 2700                                # density of the beam [kg/m^3]

# List definitions
lst_new_horimotor = [[] for x1 in range(5)]
lst_new_horiprop = [[] for x2 in range(5)]
lst_new_motor = [[] for x3 in range(5)]
lst_new_prop = [[] for x4 in range(5)]
lst_m_motor_structure = [[] for x5 in range(5)]
lst_m_battery = [[] for x6 in range(5)]

# Iteration for mass and power for all concepts
for i in range(0, 5):
    for j in range(0, 5):

        # ----HOVER POWER per 1 rotor------
        total_T = mass_f * g                                # total thrust needed to lift vehicle [N]
        Total_A = total_T / disk_loading[i]                    # total rotor area vehicle [m^2]
        one_rotor_area = Total_A / rotors_number[i]            # rotor area of ONE rotor [m^2]
        rotor_d = np.sqrt(4 * one_rotor_area / np.pi)       # rotor diameter [m]

        P_hover = FOM*np.sqrt((total_T / rotors_number[i])** 3 / (rho * one_rotor_area))  # Power required for hovering [W]
        # print("Hover Power = ", P_hover/1000, "kW")

        # ------CLIMB/DESCEND POWER required ONE rotor----
        v_i = np.sqrt((total_T / rotors_number[i]) / (2 * rho * one_rotor_area))  # ???????????????????????????????????

        P_climb = P_hover                                       # Climb power required ONE rotor [W]
        P_descend = P_hover                                     # Descend power required ONE rotor [W]

        # ------CRUISE POWER------##
        D = 1.1 * (0.5 * rho * (V_cruise) ** 2 * S * CD0)       # Drag [N] -> 1.1 for -5alpha drag
        P_cruise = P_hover + (D * V_cruise / rotors_number[i])     # Power required during cruise [W]

        # print("------------")
        # print("P_hover per 1 rotor", P_hover / 1000, "kW")
        # print("P_climb per 1 rotor", P_climb / 1000, "kW")
        # print("P_descend per 1 rotor", P_descend / 1000, "kW")
        # print("P_cruise per 1 rotor", P_cruise / 1000, "kW")
        # print("------------")

        # mass of a motor [kg]
        m_motor = (0.188 * rotors_number[i] * P_cruise / 1000 + 5.836) / rotors_number[i]  # power in the equation must be given in kW thus P_climb/1000 (article figure)

        # mass of a propeller [kg]
        m_prop = 1.1 * (rotor_d * (P_cruise * rotors_number[i] / 1000 / rotors_number[i]) * np.sqrt(blades_number)) ** 0.52

        # mass of a strut [kg]
        m_motor_structure = structure_length[i] * structure_Area * structure_density

        # mass of a battery [kg]
        m_battery = rotors_number[i] * 2 * (
                    P_hover * t_hover + P_climb * t_climb + P_cruise * t_cruise + P_descend * t_descend) / \
                    (battery_density * 3600 * battery_efficiency)
        #
        # print("motor structure mass", m_motor_structure)
        # print("rotor mass", m_motor)
        # print("propeller mass", m_prop)
        # print("battery mass", m_battery)

        # list of iterations for different mass of rotors and propellers
        lst_new_motor[i].append(m_motor)
        lst_new_prop[i].append(m_prop)
        lst_m_motor_structure[i].append(m_motor_structure)
        lst_m_battery[i].append(m_battery)

        # update aircraft j ==0 and then j> exchange the rotors and propellors
        if j == 0:
            mass_f = mass_f + rotors_number[i] * (m_motor + m_prop) + m_battery + 4 * structure_penalty[i] * m_motor_structure
        else:
            mass_f = mass_f + rotors_number[i] * (lst_new_prop[i][j] + lst_new_motor[i][j]) + 4 * structure_penalty[i] * \
                   lst_m_motor_structure[i][j] - rotors_number[i] * (
                               lst_new_prop[i][j - 1] + lst_new_motor[i][j - 1]) - 4 * structure_penalty[i] * lst_m_motor_structure[i][
                       j - 1] + lst_m_battery[i][j] - lst_m_battery[i][j - 1]

        # print("mass", mass_f)
        propeller_radius = np.sqrt(one_rotor_area / np.pi)                  # radius of the propeller blades [m]
        CD0 = drag.Cd0_design1(rotors_number[i], propeller_radius, total_T)    # zero-lift drag coefficient [-]

    print("Battery Mass concept" + str(i + 1) + " = " + str(lst_m_battery[i][-1]))
propeller_radius = np.sqrt(one_rotor_area / np.pi)
print("propeller radius: ", propeller_radius, "m")

# TOTAL ENERGY
Total_Energy = rotors_number[0] * 2 * (P_hover * t_hover + P_climb * t_climb + P_cruise * t_cruise + P_descend * t_descend)
print("Concept 1 Total Energy per mission:", Total_Energy / 1000, "KJ")

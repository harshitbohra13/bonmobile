print("hello world")
import numpy as np
#data

mass = 500 #[kg] #aircraft mass
g = 9.8 
rho = 1.225
FOM = 0.75
CD0=0.3
S=7

V_descend=7 #m/s
V_climb = 7 #m/s
V_cruise=100/3.6 #m/s

t_climb=100 #s
t_descend=100 #s
t_cruise = 20/100 *3600 + 120 #[s] +120 seconds for acceleration and 20 km in one direction 
t_hover = 60 #s

battery_efficiency = 0.85
battery_density = 250 #Wh/kg

#---concept 1---
#rotors data
blades_number = 3
rotors_number = 12
disk_loading= 500 # [N/m^2]
#circular beam
structure_penalty = 1.5 #structure penalty for additional mass with regards to concept 1,
structure_length = 2
structure_radius = 0.075
structure_inner_radius = 0.070
structure_Area = np.pi*structure_radius**2 - np.pi*structure_inner_radius**2
structure_density = 2700 #[kg/m**3]
#-------concept 2------------
#concept 3
#concept 4




#list definition
lst_new_motor = []
lst_new_prop = []
lst_m_motor_structure = []
lst_m_battery = []

#iteration for mass and power
i= 0
for i in range(0,5):
    
    #----HOVER POWER per 1 rotor------
    total_T = mass*g # total thrust
    Total_A=total_T/disk_loading
    one_rotor_area=Total_A/rotors_number
    rotor_d= np.sqrt(4*one_rotor_area/np.pi)

    P_hover = FOM * np.sqrt((total_T/rotors_number)**3/(rho*one_rotor_area)) #[W]
    # print(P_hover, "W")
    # print()
    print(P_hover/1000, "kW")
    
    # print()

    #------CLIMB/DESCEND POWER per 1 rotor----
    
    v_i = np.sqrt((total_T/rotors_number)/(2*rho*one_rotor_area))    #???????????????????????????????????
    #v_h=np.sqrt(V_climb*np.sqrt((total_T/rotors_number))+((total_T/rotors_number)/(2*rho*one_rotor_area)))
    # print("v_h", v_h)

    P_climb= P_hover #per 1 rotor
    P_descend =  P_hover #per 1 rotor

    #------CRUISE POWER------
    D=0.5*rho* (V_cruise)**2 * S * CD0
    P_cruise=P_hover+(D*V_cruise/rotors_number)

    print("------------")
    print("P_hover per 1 rotor", P_hover/1000, "kW")
    print("P_climb per 1 rotor", P_climb/1000, "kW")
    print("P_descend per 1 rotor", P_descend/1000, "kW")
    print("P_cruise per 1 rotor", P_cruise/1000, "kW")
    print("------------")



    #mass of a motor
    m_motor = (0.188*rotors_number*P_climb/1000 +5.836)/rotors_number #power in the equation must be given in kW thus P_climb/1000 (article figure)
    #mass 0f a propeller
    m_prop = 1.1 *(rotor_d*(P_climb*rotors_number/1000/rotors_number)*np.sqrt(blades_number))**0.52
    #mass of a strut
    m_motor_structure = structure_length * structure_Area * structure_density
    #mass of a battery
    m_battery = rotors_number*2*(P_hover*t_hover +  P_climb*t_climb + P_cruise*t_cruise + P_descend * t_descend)/(battery_density*3600*battery_efficiency)

    print("motor structure mass",m_motor_structure)
    print( "rotor mass", m_motor)
    print( "propeller mass", m_prop)
    print("battery mass", m_battery)

    #list of iterartions for different mass of rotors and propellers
    lst_new_motor = lst_new_motor + [m_motor]
    lst_new_prop = lst_new_prop+ [m_prop]
    lst_m_motor_structure = lst_m_motor_structure + [m_motor_structure]
    lst_m_battery = lst_m_battery + [m_battery]
    #update aircraft i ==0 and then i> exchange the rotors and propellors
    if i ==0:
        mass = mass + rotors_number*(m_motor+m_prop)+m_battery + 4*structure_penalty * m_motor_structure
    else:
        mass = mass + rotors_number*(lst_new_prop[i] + lst_new_motor[i])+4*structure_penalty *lst_m_motor_structure[i] - rotors_number*(lst_new_prop[i-1] + lst_new_motor[i-1])- 4* structure_penalty *lst_m_motor_structure[i-1] + lst_m_battery[i]-lst_m_battery[i-1]
    print("mass", mass)
    i=i+1
propeller_radius = np.sqrt(one_rotor_area/np.pi)
print("propeller radius: ", propeller_radius,"m")

#TOTAL ENERGY:
Total_Energy = rotors_number*2*(P_hover*t_hover +  P_climb*t_climb + P_cruise*t_cruise + P_descend * t_descend)
print("Concept 5 Total Energy per mission:",Total_Energy/1000, "KJ" )

#check again:
#penalty factor 10 struts, 20 engines


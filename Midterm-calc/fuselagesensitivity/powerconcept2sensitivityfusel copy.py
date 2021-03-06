print("hello world")
import numpy as np
import drag_calc as drag
import matplotlib.pyplot as plt

#airfoil charachteristics
def sensfus(massfus):  
    alpha=5
    C_D=0.008
    C_L=0.5


    S_airfoil=10
    n_ult=3.75
    AR=8.3

    # mass = 500+ 0.04674*(1000**0.397)*(S_airfoil**0.36)*(n_ult**0.397)*AR**1.712 #[kg] #aircraft 
    mass = massfus
    # Sw = 20
    # t_c = 0.15
    # taper_ratio = 0.3
    # mass_WING = (0.0051*(1200*n_ult)**0.557*Sw**0.649*AR**0.5*t_c**(-0.4)*(1+taper_ratio)**0.1*np.cos(np.radians(20))**(-1)*1)#[kg] #aircraft mass
    # print("!!!!!!!!!!!!!!",mass_WING)
    W_to = 1490 * 9.8
    W_D = W_to
    n_ult=3.75
    S = S_airfoil
    b = np.sqrt(AR*S)
    # print("wingspan:",b)
    t_c = 0.12
    taper_ratio = 0.7
    f_lambda = (1+ 2* taper_ratio)/(1+taper_ratio)
    sweep = 0
    specific_material_weight = 2700*9.8 #N/m^3
    t_ss = 0.004 * (1 + np.sqrt(W_to/(10**6)))
    w_sigma = 40 * (1+1.1*(W_to/(10**6))**(-1/4))*10**(-6)
    Weight_WING = 0.06 * w_sigma * n_ult * W_D * b**3 / S * f_lambda / (t_c * np.cos(np.radians(sweep))**2) + specific_material_weight*t_ss*S

    g = 9.8 
    mass = mass+ Weight_WING/g #[kg] #aircraft mass
    # print(mass-500)
    print("WIIIIIIING:", Weight_WING/g)

    rho = 1.225
    FOM = 0.75
    CD0=0.3
    S=7 #???????
    # Weight_WING = 100
    # mass = mass+ Weight_WING 


    V_descend=7 #m/s
    V_climb = 7 #m/s
    V_cruise=100/3.6 #m/s

    t_climb=(450-30.5)/ V_climb #s
    t_descend= t_climb #s
    t_cruise = 20/100 *3600 + 120 #[s] +120 seconds for acceleration and 20 km in one direction 
    t_hover = 60 *2 #[s] hovering appears twice

    #BATTERIES
    battery_efficiency = 0.85
    battery_density = 250 #Wh/kg
    #HYDROGEN
    # battery_efficiency = 0.9
    # battery_density = 33000 #Wh/kg

    #---concept 1---
    #rotors data
    blades_number = 3
    rotors_number = 4
    disk_loading= 500 # [N/m^2]
    #circular beam
    structure_penalty = 1 #structure penalty for additional mass with regards to concept 1
    structure_length = 1
    structure_radius = 0.075
    structure_inner_radius = 0.070
    structure_Area = np.pi*structure_radius**2 - np.pi*structure_inner_radius**2
    structure_density = 2700 #[kg/m**3]
    #-------concept 2------------
    #concept 3
    #concept 4




    #list definition
    lst_new_horimotor = []
    lst_new_horiprop = []
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
        P_hovercruise=FOM * np.sqrt(((total_T-(0.5*rho*V_cruise**2*C_L*S_airfoil))/rotors_number)**3/(rho*one_rotor_area))
        #------CRUISE POWER------
        D=0.5*rho* (V_cruise)**2 * S * CD0
        D_airfoil=0.5*rho*(V_cruise)**2*S_airfoil*C_D 
        Total_drag=D_airfoil+D

        P_totalcruise=P_hovercruise+Total_drag*V_cruise

        print("------------")
        print("P_hover per 1 rotor", P_hover/1000, "kW")
        print("P_climb per 1 rotor", P_climb/1000, "kW")
        print("P_descend per 1 rotor", P_descend/1000, "kW")
        print("P_cruise per 1 rotor", P_hovercruise/1000, "kW")

        # print("P_cruise per 1 rotor thrust",(total_T-(0.5*rho*V_cruise**2*C_L*S_airfoil))/rotors_number, "N")
        print("------------")

        # horizontal propellor calc
        horprop_area=Total_drag/disk_loading
        horprop_d= np.sqrt(4*horprop_area/np.pi)
        m_horiprop=1.1 *(horprop_d*((Total_drag*V_cruise)/1000)*np.sqrt(blades_number))**0.52
        m_horimotor = (0.188*(Total_drag*V_cruise)/1000 +5.836)

        #mass of a motor
        m_motor = (0.188*rotors_number*P_climb/1000 +5.836)/rotors_number #power in the equation must be given in kW thus P_climb/1000 (article figure)
        #mass 0f a propeller
        m_prop = 1.1 *(rotor_d*(P_climb*rotors_number/1000/rotors_number)*np.sqrt(blades_number))**0.52
        #mass of a strut
        m_motor_structure = structure_length * structure_Area * structure_density
        #mass of a battery
        m_battery = (2*Total_drag*V_cruise*t_cruise + rotors_number*2*(P_hover*t_hover +  P_climb*t_climb + P_hovercruise*t_cruise + P_descend * t_descend))/(battery_density*3600*battery_efficiency)

        print("motor structure mass",4*m_motor_structure)
        print( "rotor mass", m_motor)
        print( "propeller mass", m_prop)
        print( "hor_rotor mass", m_horimotor)
        print( "hor propeller mass", m_horiprop)
        print("battery mass", m_battery)
        print("hor propeller radius", horprop_d/2)

        #list of iterartions for different mass of rotors and propellers
        lst_new_motor = lst_new_motor + [m_motor]
        lst_new_prop = lst_new_prop+ [m_prop]
        lst_m_motor_structure = lst_m_motor_structure + [m_motor_structure]
        lst_m_battery = lst_m_battery + [m_battery]
        lst_new_horimotor =  lst_new_horimotor +  [m_horimotor]
        lst_new_horiprop = lst_new_horiprop + [m_horiprop]
        #update aircraft i ==0 and then i> exchange the rotors and propellors
        if i ==0:
            mass = mass + rotors_number*(m_motor+m_prop)+m_battery + 4*structure_penalty * m_motor_structure + (m_horiprop+m_horimotor)
        else:
            mass = mass + rotors_number*(lst_new_prop[i] + lst_new_motor[i])+4*structure_penalty *lst_m_motor_structure[i] + (lst_new_horimotor[i]+lst_new_horiprop[i]) - rotors_number*(lst_new_prop[i-1] + lst_new_motor[i-1])- 4* structure_penalty *lst_m_motor_structure[i-1] - (lst_new_horimotor[i-1]+lst_new_horiprop[i-1]) + lst_m_battery[i]-lst_m_battery[i-1]
        print("mass", mass)
        
        propeller_radius = np.sqrt(one_rotor_area/np.pi)
        CD0 = drag.Cd0_design2(rotors_number, propeller_radius, total_T, P_hovercruise/V_cruise, horprop_d)
        i=i+1
    propeller_radius = np.sqrt(one_rotor_area/np.pi)
    print("propeller radius: ", propeller_radius,"m")

    #TOTAL ENERGY:
    Total_Energy = (2*Total_drag*V_cruise*t_cruise + rotors_number*2*(P_hover*t_hover +  P_climb*t_climb + P_hovercruise*t_cruise + P_descend * t_descend))
    print("Concept 2 Total Energy per mission:",Total_Energy/1000, "KJ" )


    # print()
    # print("Concept 2 Total Energy per mission:",Total_Energy/10**6, "MJ" )
    # print("Number of rotors:", " 4 rotors with vertical thrust, 1 rotor with horizontal thrust")
    # print("------------")
    # print("Mass of each vertical thrust rotor", m_motor)
    # print("Mass of each horizontal thrust rotor", m_horimotor)
    # print( "propeller mass for vertical thrust rotor", m_prop)
    # print("propeller radius vertical thrust: ", propeller_radius,"m")
    # print( "propeller mass for horizontal thrust rotor", m_horiprop)
    # print("propeller radius horizontal thrust: ", horprop_d/2,"m")
    # print("------------")
    # print("Wing mass: ", Weight_WING/g)
    # print("motor structure mass",4*m_motor_structure)
    # print("battery mass", m_battery)
    # print("Total mass", mass)
    # print("------------")
    # print("P_hover per 1 rotor vertical thrust", P_hover/1000, "kW")
    # print("P_climb per 1 rotor vertical thrust", P_climb/1000, "kW")
    # print("P_descend per 1 rotor vertical thrust", P_descend/1000, "kW")
    # print("P_cruise per 1 rotor vertical thrust", P_hovercruise/1000, "kW")
    # print("P_cruise per 1 rotor horizontal thrust", Total_drag*V_cruise/1000, "kW")
    # print("------------")
    return mass


lstfus=np.arange(280,501,1)
lst_totalmass=[]
for i in lstfus:
    lst_totalmass.append(sensfus(i))


print(lstfus)
print(lst_totalmass)
plt.plot(lstfus,lst_totalmass)
plt.xlabel("Fuselage mass [kg]")
plt.ylabel("Total mass [kg]")
plt.show()
from re import A
import numpy as np
#data

mass = 700 #[kg] #aircraft mass
g = 9.8 
rho = 1.225
FOM = 0.75

V_climb = 5
v_h = 1

#rotors data
#concept 1
rotor_d = 1.5 #[m]
blades_number = 3
rotors_number = 4
Area = np.pi *(rotor_d/2)**2 #[m^2] #propeler area
#concept 2
#concept 3
#concept 4




#list definition
lst_new_motor = []
lst_new_prop = []

#iteration for mass and power
i= 0
for i in range(0,5):
    
    #HOVER POWER
    T = mass*g #thrust
    P_hover = FOM * np.sqrt(T**3/(rho*Area)) #[W]
    # print(P_hover, "W")
    print()
    print(P_hover/1000, "kW")
    print()

    #CLIMB POWER
    # P_climb = P_hover
    P_climb = P_hover*(V_climb/(2*v_h) + np.sqrt((V_climb/(2*v_h))**2+1))
    P_descend =  P_hover*(V_climb/(2*v_h) - np.sqrt((V_climb/(2*v_h))**2-1))

    #mass of a motor
    m_motor = (0.188*P_climb/1000 +5.836)/rotors_number #power in the equation must be given in kW thus P_climb/1000 (article figure)
    #mass 0f a propeller
    m_prop = 1.1 *(rotor_d*(P_climb/1000/rotors_number)*np.sqrt(blades_number))**0.52

    print( "rotor mass", m_motor)
    print( "propeller mass", m_prop)

    #list of iterartions for different mass of rotors and propellers
    lst_new_motor = lst_new_motor + [m_motor]
    lst_new_prop = lst_new_prop+ [m_prop]

    
    #update aircraft i ==0 and then i> exchange the rotors and propellors
    if i ==0:
        mass = mass + rotors_number*(m_motor+m_prop)
    else:
        mass = mass + rotors_number*(lst_new_prop[i] + lst_new_motor[i]) - rotors_number*(lst_new_prop[i-1] + lst_new_motor[i-1]) 
    print("mass", mass)
    i=i+1




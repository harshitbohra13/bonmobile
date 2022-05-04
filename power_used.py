#times for particular phases
# from matplotlib.font_manager import _Weight

import numpy as np
# t_climb =
# t_descend = 
# t_cruise = 40/100 #range over crusie speed

#data
mass = 700 #[kg] #aircraft mass
g = 9.8 
rho = 1.225
rotor_d = 1.5 #[m]

#rotors data
blades_number = 3
rotors_number = 4
Area = np.pi *(rotor_d/2)**2 #[m^2] #propeler area
FOM = 0.75

# T = mass*g #thrust [N]

#iteration 
i= 0
for i in range(0,5):
    i=i+1
    T = mass*g
    P_hover = FOM * np.sqrt(T**3/(rho*Area)) #[W]
    # print(P_hover, "W")
    print(P_hover/1000, "kW")

    P_climb = P_hover

    m_motor = (0.188*P_climb/1000 +5.836)/rotors_number #power must be given in kW article figure
    m_prop = 1.1 *(rotor_d*(P_climb/1000/rotors_number)*np.sqrt(blades_number))**0.52
    print( "rotor mass", m_motor)
    print( "propeller mass", m_prop)
    
    m_new_motor = m_motor
    m_new_prop = m_prop

    mass = mass + rotors_number*(m_motor+m_prop)
    print("mass", mass)




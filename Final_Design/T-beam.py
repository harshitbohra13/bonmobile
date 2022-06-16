import numpy as np
Force_Z = 7185 #N 7185 
Moment_y = 62.416 #62.416 
Moment_x = 1191.466 #1191.466

#iteration 0 
Force_Z = 7263 #N
Moment_y = 28.238 #62.416 
Moment_x = 1235.8 #1191.466

# #iteration 1
# Force_Z = 8729 #N
# Moment_y = 33.9 #62.416 
# Moment_x = 1485.4 #1191.466


Density  = 2770
yielding_strength = 469 #Mpa
safety_factor = 1.5
yielding_strength_with_safety_fact = 469/safety_factor

B_f = 47 #mm
T_f = 3 #mm
H = 50 #mm
T_w = 3 #mm

#89%


#Centeroid zero at T_W level where H starts

A_1 = T_w*H
A_1_center = H/2
A_2 = B_f*T_f
A_2_center = (H+T_f/2)

#Centroid
Centroid = (A_1*A_1_center+A_2*A_2_center)/(A_1+A_2)/1000 #m
print(Centroid,"m")
#Sec moment of area

Force_Z_stress = Force_Z / (A_1+A_2)
print("Stress due to Z force",Force_Z_stress, "MPa")

#XX moment of area
I_xx  = 1/12 * T_w*H**3 + A_1*(A_1_center-Centroid*1000)**2 + 1/12*B_f*T_f**3+A_2*(A_2_center-Centroid*1000)**2
I_xx = I_xx/(10**12)
# I_xx =213680/(10**12)
#YY moment of area
I_yy = 1/12 * H * T_w**3 + 1/12* T_f*B_f**3
I_yy = I_yy/(10**12)
print(I_xx,"Ixx")
print(I_yy,"Iyy")

#stresses 
Moment_x_stress = abs( Moment_x * (((H+T_f)/1000)-Centroid)/I_xx)/(10**6)
Moment_y_stress = abs(Moment_y* ((B_f/2)/1000)/I_yy)/(10**6)
print("Moment x stress",Moment_x_stress)
print("Moment y stress",Moment_y_stress)
Max_stress = Force_Z_stress+Moment_x_stress+Moment_y_stress

print()
print("Max_stress",Max_stress,"MPa")
print("Maximal allowable stress with SF",yielding_strength_with_safety_fact,"MPa")
print()
Circ1 = np.pi*(2.082/2+1.550/2)
Circ2 = np.pi*(2.550/2+1.900/2)
Circ3 = np.pi*(2.550/2+1.900/2)
Circ4 = np.pi*(2.410/2+1.800/2)

print("Circ_1:",Circ1)
print("Circ_2:",Circ2)
print("Circ_3:",Circ3)
print("Circ_4:",Circ4)
Total_Volume = (Circ1+Circ2+Circ3+Circ4)*(A_1+A_2)/10**6
print("per:",Circ1+Circ2+Circ3+Circ4)
ribs_mass = Total_Volume * Density
print("Mass",ribs_mass)
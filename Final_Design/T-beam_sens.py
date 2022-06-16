from cProfile import label
import numpy as np
import matplotlib.pyplot as plt

sigma_cr_lst = []
H_lst = []
B_f = 47 #mm
T_f = 3 #mm
H = 25 #mm
T_w = 3 #mm

for i in range(0,50):
    #iteration 0 
    Force_Z = 7263 #N
    Moment_y = 28.238 #62.416 
    Moment_x = 1235.8 #1191.466


    Density  = 2770
    yielding_strength = 469 #Mpa
    safety_factor = 1.5
    yielding_strength_with_safety_fact = 469/safety_factor



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

    sigma_cr_lst = sigma_cr_lst + [Max_stress]
    H_lst = H_lst + [H]
    H = H+1

print(sigma_cr_lst)

#graph
# x axis values
x = H_lst
# corresponding y axis values
y = sigma_cr_lst
  
# plotting the points 
plt.plot(x, y, label = "Critical stress in the cross-section")
plt.plot([25,75], [312,312], label="Allowable stress")
# naming the x axis
plt.xlabel('H dimension [mm]')
# naming the y axis
plt.ylabel('Critical stress [MPa]')
  
# giving a title to my graph
#plt.title('Sensitivity analysis')
plt.grid( linestyle = '--', linewidth = 0.5)
# function to show the plot
plt.legend(loc="upper right")
plt.savefig('H_dimension_sens.png')
plt.show()


#-----------------------------BF------------------------------------
sigma_cr_lst = []
B_F_lst = []
B_f = 22 #mm 47
T_f = 3 #mm
H = 50 #mm
T_w = 3 #mm

for i in range(0,50):
    #iteration 0 
    Force_Z = 7263 #N
    Moment_y = 28.238 #62.416 
    Moment_x = 1235.8 #1191.466


    Density  = 2770
    yielding_strength = 469 #Mpa
    safety_factor = 1.5
    yielding_strength_with_safety_fact = 469/safety_factor



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

    sigma_cr_lst = sigma_cr_lst + [Max_stress]
    B_F_lst = B_F_lst + [B_f]
    B_f = B_f+1

print(sigma_cr_lst)

#graph
# x axis values
x = B_F_lst
# corresponding y axis values
y = sigma_cr_lst
  
# plotting the points 
plt.plot(x, y, label = "Critical stress in the cross-section")
plt.plot([22,72], [312,312], label="Allowable stress")
# naming the x axis
plt.xlabel('B dimension [mm]')
# naming the y axis
plt.ylabel('Critical stress [MPa]')
  
# giving a title to my graph
#plt.title('Sensitivity analysis')
plt.grid( linestyle = '--', linewidth = 0.5)
# function to show the plot
plt.legend(loc="upper right")
plt.savefig('B_dimension_sens.png')
plt.show()

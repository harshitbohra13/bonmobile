import numpy as np
# The design is based on the 3D diagam and hollow cylinders as rods

#---------------dimensions-----------------
#configuration/ rod geometry
h = 0 # vertical height of the rod
x = 0 # horizontal distance of the rod end to fuselage
r_out_rod = 0.05
r_in_rod = 0.03
#ski
ski_width = 0
ski_length = 0
ski_thickness = 0
d_inner = 0
d_outer = 0

#-------------material----------------
#rod
mat_dens_rod = 0
mat_strength_rod = 0
#ski
mat_dens_ski= 0
mat_strength_ski = 0

safety_factor = 1.5
#-----------force in ONE rod--------------
Max_Weight = 0
F_in_rod = np.sin(np.arctan(h/x)) * Max_Weight/4


#----------changing the inner radius to fit the rod stress-----------------
i = True
while i == True:

#Cross-section
Area_skid_rod = np.pi*(r_out_rod**2-r_in_rod**2)
#stress is euqally distributed over the cross-section
stress_rod = F_in_rod / Area_skid_rod
if stress_rod > mat_strength_rod / safety_factor:
    r_in_rod = r_in_rod - r_in_rod - 0.001
else:
    i = False

#Moment of intertia cross-section
I_yy = math.pi/64* (d_outer^4-d_inner^4) 
#----------------calculating mass-------------------------
rod_mass = 4 * Area_skid_rod * np.sqrt(h**2+x**2) * mat_dens_rod
ski_mass = 2 * ski_width * ski_length * ski_thickness * mat_dens_ski

skid_mass = rod_mass+ski_mass
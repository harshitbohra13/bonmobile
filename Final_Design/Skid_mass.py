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

#----------------calculating mass-------------------------
rod_mass = 4 * Area_skid_rod * np.sqrt(h**2+x**2) * mat_dens_rod
ski_mass = 2 * ski_width * ski_length * ski_thickness * mat_dens_ski

skid_mass = rod_mass+ski_mass

#-----------forces and moments at the root --------------
F_x_end = 0
F_z_end = 0
M_x_root = F_z_end * x
M_y_root = F_x_end * h
M_z_root = F_x_end * x

M_z_1 = M_z * sin(alpha)
M_z_2 = M_z * cos(alpha)
M_y_1 = M_y * cos(alpha)
M_y_2 = M_y * sin(alpha)



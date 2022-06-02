import numpy as np
# The design is based on the 3D diagam and hollow cylinders as rods
h = 0 # vertical height of the rod
x = 0 # horizontal distance of the rod end to fuselage
ski_width = 0
ski_length = 0
ski_thickness = 0
mat_dens = 0

#force in the rod
Max_Weight = 0
F_in_rod = np.sin(np.arctan(h/x)) * Max_Weight/4

#rod geometry
r_out_rod = 0.05
r_in_rod = 0.03

Area_skid_rod = np.pi*(r_out_rod**2-r_in_rod**2)

#stress is euqally distributed over the cross-section
stress_rod = F_in_rod / Area_skid_rod

mat_dens = 
# rod_mass = 

#find mass = 
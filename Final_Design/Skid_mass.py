import numpy as np
# The design is based on the 3D diagam and hollow cylinders as rods

#------- INPUT BY USER // test parameters -------
#configuration/ rod geometry
h = 0 # vertical height of the rod [m]
w = 0 # width rod [m]
r_in_rod = 0.03 # inner diamater rod [m]
F_x_rod = 1000 #[N]
F_z_rod = 10000 #[N]
#---------------dimensions-----------------
#configuration/ rod geometry
# h = h # vertical height of the rod [m]
#w = w # width rod [m]

l = np.sqrt(h*h + w*w) # length rod
alpha = np.arctan(w/h) # angle rod

r_out_rod = r_in_rod + t # outer diameter rod [m]
#r_in_rod = 0.03 # inner diamater rod [m]

#ski
ski_width = 0
ski_length = 0
ski_thickness = 0


#-------------material----------------
#rod
mat_dens_rod = 2710 #kg/m^3
mat_strength_rod = 276 #MPa
mat_yield_rod = 276 #MPa
mat_E_modulus = 72.4 #GPa
mat_G_modulus = 24 #GPa

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
mass_rod = mat_dens_rod * np.pi * (r_out_rod**2 - r_in_rod**2) * h/cos(alpha)   #mass rod
W_rod = mass_rod * 9.81 # weight rod

R_x = F_x_end
R_y = 0
R_z = W_rod - F_z_end
M_x = W_rod * 0.5 * w - F_z_end * w
M_y = -F_x_end * h
M_z = -F_x_end * w

#----------stresses in cross-section ----------
t = 0
while t < r_in_rod:
    sigma_normal = R_z / (np.pi * t**2 + 2 * t * r_in_rod)
    sigma_shear = R_x * (t + r_in_rod) * 0.5 * np.pi * (t**2 + 2 * t * r_in_rod) / (np.pi/4 * (t**4 + 4 * t**3 * r_in_rod + 6*t**2 * r_in_rod**2 + 4*t*r_in_rod**3) * t)
    sigma_torsion = M_z * (t + r_in_rod) / ( np.pi/2 * (t**4 + 4 * t**3 * r_in_rod + 6*t**2 * r_in_rod**2 + 4*t*r_in_rod**3))
    sigma_bending = (M_x + M_y) * (t + r_in_rod) / (np.pi/4 * (t**4 + 4 * t**3 * r_in_rod + 6*t**2 * r_in_rod**2 + 4*t*r_in_rod**3))

    stresses  = sigma_normal + sigma_shear + sigma_torsion + sigma_bending
    max_stress = mat_yield_rod *100000 #[Pa]

    if stresses > max_stress:
        print("thickness rod is", t, ", and total stress in cross section is", stresses, "Pa")
        break

    t = t+0.001

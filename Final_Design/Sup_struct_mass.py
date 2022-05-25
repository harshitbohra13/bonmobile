from sqlalchemy import null
from Battery_mass import Max_Total_thrust
from configuration_and_materials import Rotor_number

# Total_Supporting_structure_mass = 0


#Thrust per rotor
Max_thrust_per_rotor = Max_Total_thrust/Rotor_number

def Strut_1(Fuselage_width, fuselage_rotor_clearance, rotor_radius, safety_factor, T_max, r_outer, Density_1):
    Length_1 = (Fuselage_width/2) * np.sqrt(2) + np.sqrt(2) * (fuselage_rotor_clearance+rotor_radius)
    max_stress = y_strength_mat1/safety_factor
    M = 0.25* T_max * Length_1
    I = (M*r_outer)/(max_stress)
    r_inner = (r_outer**4 - (4/np.pi)*I)**0.25
    Area_1 = np.pi *(r_outer**2 - r_inner **2)
    Strut_1_mass = Area_1 * Length_1 * Density_1
    return Strut_1_mass
    # print("max stress", max_stress/(10**6), "MPa")
    
def sup_struct_mass(Rotor_number, ):
    return null
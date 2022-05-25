from sqlalchemy import null
from Battery_mass import Max_Total_thrust
from configuration_and_materials import configure, Aluminiumbadboy
from Fuselage_mass import *
from Motor_mass import *
import numpy as np

# Total_Supporting_structure_mass = 0


#Thrust per rotor
def thrustperrotor(Total_T, Rotor_number):
    Max_thrust_per_rotor = Max_Total_thrust/Rotor_number
    return(Max_thrust_per_rotor)

def Strut_1(bon, Fuselage_width, T_max, r_outer):
    Length_1 = (Fuselage_width/2) * np.sqrt(2) + np.sqrt(2) * (bon.fuselage_rotor_clearance+bon.rotor_radius)
    max_stress = Aluminiumbadboy.y_strength_mat/bon.safety_factor
    M = 0.25* T_max * Length_1 + 2*T_max* bon.rotor_to_rotor_clearance 
    I = (M*r_outer)/(max_stress)
    r_inner = (r_outer**4 - (4/np.pi)*I)**0.25
    Area_1 = np.pi *(r_outer**2 - r_inner **2)
    Strut_1_mass = Area_1 * Length_1 * Aluminiumbadboy.Density
    return Strut_1_mass
    # print("max stress", max_stress/(10**6), "MPa")
    
def sup_struct_mass(Rotor_number, ):
    return null

def Strut_2(bon, T_max, r_outer):
    Length_2 = 2* bon.rotor_radius + bon.rotor_to_rotor_clearance
    max_stress = Aluminiumbadboy.y_strength_mat/bon.safety_factor
    M =  (T_max * Length_2)/12
    I = (M*r_outer)/(max_stress)
    r_inner = (r_outer**4 - (4/np.pi)*I)**0.25
    Area_2 = np.pi *(r_outer**2 - r_inner **2)
    Strut_2_mass = Area_2 * Length_2 * Aluminiumbadboy.Density
    return Strut_2_mass
    # print("max stress", max_stress/(10**6), "MPa")

def Strut_3(bon, T_max, r_outer):
    Length_2 = 1
    max_stress = Aluminiumbadboy.y_strength_mat/bon.safety_factor
    M =  1
    I = (M*r_outer)/(max_stress)
    r_inner = (r_outer**4 - (4/np.pi)*I)**0.25
    Area_2 = np.pi *(r_outer**2 - r_inner **2)
    Strut_2_mass = Area_2 * Length_2 * Aluminiumbadboy.Density
    return Strut_2_mass
    # print("max stress", max_stress/(10**6), "MPa")
# from sqlalchemy import null
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

def Strut_1(bon):
    Length_1 = (bon.fuselage_width/2) * np.sqrt(2) + np.sqrt(2) * (bon.fuselage_rotor_clearance+bon.rotor_radius)
    max_stress = Aluminiumbadboy.y_strength_mat/bon.safety_factor
    M = 0.25 * bon.T_max * Length_1 + 2*bon.T_max* bon.rotor_to_rotor_clearance
    I = (M*bon.r_outer)/(max_stress)
    r_inner = (bon.r_outer**4 - (4/np.pi)*I)**0.25
    Area_1 = np.pi *(bon.r_outer**2 - r_inner **2)
    Strut_1_mass = Area_1 * Length_1 * Aluminiumbadboy.Density
    return Strut_1_mass, Length_1
    # print("max stress", max_stress/(10**6), "MPa")
    
def sup_struct_mass(Rotor_number, ):
    return null

def Strut_2(bon):
    Length_2 = 2* bon.rotor_radius + bon.rotor_to_rotor_clearance
    max_stress = Aluminiumbadboy.y_strength_mat/bon.safety_factor
    M = (bon.T_max * Length_2)/12
    I = (M*bon.r_outer)/(max_stress)
    r_inner = (bon.r_outer**4 - (4/np.pi)*I)**0.25
    Area_2 = np.pi *(bon.r_outer**2 - r_inner **2)
    Strut_2_mass = Area_2 * Length_2 * Aluminiumbadboy.Density
    return Strut_2_mass, Length_2
    # print("max stress", max_stress/(10**6), "MPa")

def Strut_3(bon):
    Length_2 = 2* bon.rotor_radius + bon.rotor_to_rotor_clearance
    max_stress = Aluminiumbadboy.y_strength_mat/bon.safety_factor
    M = (bon.T_max * Length_2)/12
    I = (M*bon.r_outer)/(max_stress)
    r_inner = (bon.r_outer**4 - (4/np.pi)*I)**0.25
    Area_2 = np.pi *(bon.r_outer**2 - r_inner **2)
    Strut_2_mass = Area_2 * Length_2 * Aluminiumbadboy.Density
    return Strut_2_mass, Length_2
    # print("max stress", max_stress/(10**6), "MPa")

bon = configure("bon")
strut1 = Strut_1(bon)
strut2 = Strut_2(bon)
strut3 = Strut_3(bon)

print(strut1[1], strut2[1], strut3[1])